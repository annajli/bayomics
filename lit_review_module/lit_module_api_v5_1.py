import pandas as pd
import numpy as np

import os
import json
import time
import collections
import requests
import xml.etree.ElementTree as ET

"""
PubMed E-utils (esearch/efetch) + INDRA API wrapper functions.

--- Why this replaces pubtator_api_module_v3 ---

PubTator3's entity system only normalizes six/seven bioconcept types (gene,
disease, chemical, species, genetic variant, cell line, sometimes chromosome).
Over half of this dataset's ~531 features -- cell-subset names, Hallmark/CyTOF
pathway module scores, composite readouts, demographic covariates -- have no
PubTator entity to resolve to at all, so gating the pipeline on entity
resolution silently dropped the majority of edges. This module instead:

  1. Looks up each feature's plain-English description in the curated feature
     glossary (feature_glossary.xlsx) -- a local dictionary lookup, no API call.
  2. Uses that description as free text against PubMed's esearch/efetch
     endpoints (no entity normalization required).
  3. Relies on INDRA -- not PubTator -- for structured/typed relationship
     signal (see below). PubTator was tried as an intermediate step (bonus
     BioRED relation extraction, PMID-keyed, no entity resolution required)
     but was removed: since PubTator only grounds literal gene/chemical/
     disease terms, and only a minority of this feature set literally is one,
     it contributed near-empty signal on most edges while adding a second
     external dependency and an extra failure surface -- net noise rather
     than net evidence, once INDRA was running unconditionally on every edge.

--- PubMed E-utils ---

1) esearch: free-text query -> list of PMIDs
   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
     ?db=pubmed&term=QUERY&retmode=json&retmax=N

2) efetch: PMIDs -> full record XML (title + structured abstract)
   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
     ?db=pubmed&id=PMID1,PMID2,...&rettype=abstract&retmode=xml

NCBI asks for a tool name / email on every call and grants a higher rate limit
(10 req/s vs. 3 req/s) if you supply an API key
(https://www.ncbi.nlm.nih.gov/account/settings/). Set NCBI_API_KEY / NCBI_EMAIL
env vars -- both optional, but worth it once you're running this at full
edge-list scale (esearch + efetch = 2 calls/edge minimum).

--- INDRA ---

INDRA takes free-text subject/object/agent names directly and grounds them
internally, so it never depended on PubTator entity resolution. This is now
the pipeline's only structured-relationship signal (see get_indra_statements()
/ get_directional_indra_support() below), checked THREE ways per edge:

  - subject=from_feature, object=to_feature   (forward-directional)
  - subject=to_feature, object=from_feature   (reverse-directional)
  - agent0=from_feature, agent1=to_feature    (undirected)

The subject/object endpoint only matches INDRA statement types that HAVE a
directional role (Activation, Inhibition, Phosphorylation, IncreaseAmount,
DecreaseAmount, etc.) -- per the indra_db_service README
(https://github.com/gyorilab/indra_db/blob/master/indra_db_service/README.md),
statement types without subject/object roles (Complex, Translocation,
ActiveForm, ...) are invisible to that query shape and require the `agent`
parameter instead: "This parameter is used if the specific role of the agent
... is irrelevant, or the distinction doesn't apply to the type of Statement
of interest." Without the undirected query, a real but non-directional INDRA
relationship (e.g. two proteins forming a complex) would silently count as
"no evidence" -- checking forward/reverse counts alone undercounts what
INDRA actually knows about a pair.

    GET https://db.indra.bio/statements/from_agents?subject=X&object=Y
    GET https://db.indra.bio/statements/from_agents?agent0=X&agent1=Y

Per direct guidance from the INDRA team, an API key is NOT required for this
endpoint for our use -- despite the indra_db_service README's own body text
saying one is required elsewhere in the same document. It's still accepted as
an optional `api_key` query param (only included if INDRA_API_KEY is set), in
case the team issues one later for elevated quota or access.
"""

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
INDRA_STATEMENTS_URL = "https://db.indra.bio/statements/from_agents"

# Optional -- raises the E-utils rate limit from 3 req/s to 10 req/s and is
# requested by NCBI's usage guidelines. Neither is required for this to work.
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")


def _ncbi_params(extra):
    params = dict(extra)
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL
    params.setdefault("tool", "lit_review_capstone")
    return params


def _rate_limit_sleep():
    # Stay under NCBI's unauthenticated 3 req/s cap. Skipped once an API key
    # is set, since the higher 10 req/s limit gives more headroom.
    if not NCBI_API_KEY:
        time.sleep(0.34)


# --- feature glossary (local lookup, replaces get_entity_ID as the required step) ---

def load_feature_glossary(path=None, sheet_name="Feature Glossary"):
    """
    Load the curated feature glossary into a {feature_name: description} dict.
    Call this once at orchestrator startup and pass the result into
    get_feature_concept() -- this is a local file read, not an API call.

    Parameters
    ----------
    path : str, optional
        Path to feature_glossary.xlsx. Defaults to the FEATURE_GLOSSARY_PATH
        env var, or "./feature_glossary.xlsx" if that's unset.
    sheet_name : str
        Sheet holding the "Feature Name" / "Plain-English Description" columns.

    Returns
    -------
    dict[str, str]
        feature_name -> plain-English description.
    """
    if path is None:
        path = os.environ.get("FEATURE_GLOSSARY_PATH", "./feature_glossary.xlsx")

    df = pd.read_excel(path, sheet_name=sheet_name)
    return dict(zip(df["Feature Name"], df["Plain-English Description"]))


def get_feature_concept(feature_name, glossary):
    """
    Look up a feature's plain-English concept phrase for use as a free-text
    PubMed query -- replaces the old PubTator get_entity_ID() call. Pure local
    dict lookup, no network round trip.

    Parameters
    ----------
    feature_name : str
        Raw column name from the bnlearn edge (e.g. "pb_cd14_cdc2_tnfa_nfkb").
    glossary : dict[str, str]
        Output of load_feature_glossary().

    Returns
    -------
    str or None
        The plain-English description to use as query text, or None if the
        feature isn't in the glossary.
    """
    return glossary.get(feature_name)


def get_with_retry(url, params, retries=3, backoff=2, timeout=15):
    """
    GET with retry/backoff on transient connection drops. Returns None on
    total failure instead of re-raising the final exception, so a single
    flaky request can't crash the whole orchestrator run. Callers should
    treat a None response the same as a failed/empty response.
    """
    for attempt in range(retries):
        try:
            return requests.get(url, params=params, timeout=timeout)
        except requests.exceptions.ConnectionError:
            if attempt == retries - 1:
                print(f"get_with_retry: giving up on {url} after {retries} attempts")
                return None
            time.sleep(backoff * (attempt + 1))
    return None


def build_search_term(concept_from, concept_to):
    """
    Build the exact joint esearch query term for a pair of glossary concepts.
    Factored out of get_relevant_search_results() so the identical string can
    be stamped onto results afterward for provenance (see
    attach_query_provenance()) without redefining the query logic twice.

    Parameters
    ----------
    concept_from, concept_to : str

    Returns
    -------
    str
    """
    return f"({concept_from}) AND ({concept_to})"


def get_relevant_search_results(concept_from, concept_to, retmax=20):
    """
    Query PubMed (esearch) for publications discussing both concept phrases as
    free text -- replaces PubTator3's /search/ endpoint and its accession-ID
    requirement entirely.

    Note: this queries the two concepts JOINTLY (single AND query) rather than
    separately-then-merged -- decided in favor of precision over recall, since
    both concept phrases already come from curated glossary descriptions
    rather than raw keywords, and neither this query nor the separate/merge
    alternative can encode edge direction (PubMed's search syntax has no
    "X precedes Y" operator either way) -- direction is judged downstream by
    INDRA and the LLM prompt instead.

    Parameters
    ----------
    concept_from, concept_to : str
        Plain-English concept phrases (from get_feature_concept), e.g.
        "TNF-alpha/NF-kB signaling within CD14 conventional dendritic cell 2".
    retmax : int
        Max PMIDs to retrieve.

    Returns
    -------
    list[str]
        PMIDs matching the joint query. Empty list if either concept is
        missing or the request fails.
    """
    if not concept_from or not concept_to:
        return []

    term = build_search_term(concept_from, concept_to)
    params = _ncbi_params({"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax})

    response = get_with_retry(ESEARCH_URL, params=params)
    _rate_limit_sleep()

    if response is None or response.status_code != 200 or not response.text.strip():
        return []

    try:
        payload = response.json()
    except json.JSONDecodeError:
        return []

    return payload.get("esearchresult", {}).get("idlist", [])


def get_abstract_text(pmids, batch_size=100):
    """
    Fetch title + full abstract text for a list of PMIDs via PubMed's efetch
    endpoint -- replaces PubTator3's BioC-JSON export as the TEXT source.

    Parameters
    ----------
    pmids : list[str] or str
        One PMID or a list of PMIDs (e.g. output of get_relevant_search_results).
    batch_size : int
        PMIDs per efetch call.

    Returns
    -------
    list[dict]
        One dict per publication: {"pmid": ..., "title": ..., "abstract": ...}.
        "abstract" is "" if no AbstractText element was found. Structured
        abstracts (Background/Methods/Results/Conclusions) are joined with
        their labels into a single string.
    """
    if not pmids:
        return []
    if not isinstance(pmids, (list, tuple, set)):
        pmids = [pmids]
    pmids = [str(p) for p in pmids]

    records = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        params = _ncbi_params(
            {"db": "pubmed", "id": ",".join(batch), "rettype": "abstract", "retmode": "xml"}
        )
        response = get_with_retry(EFETCH_URL, params=params, timeout=30)
        _rate_limit_sleep()

        if response is None or response.status_code != 200 or not response.text.strip():
            continue

        try:
            root = ET.fromstring(response.text)
        except ET.ParseError:
            continue

        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//PMID")
            pmid = pmid_el.text if pmid_el is not None else None

            title_el = article.find(".//ArticleTitle")
            title = "".join(title_el.itertext()).strip() if title_el is not None else ""

            abstract_parts = []
            for abst in article.findall(".//Abstract/AbstractText"):
                label = abst.get("Label")
                text = "".join(abst.itertext()).strip()
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
            abstract = " ".join(part for part in abstract_parts if part).strip()

            records.append({"pmid": pmid, "title": title, "abstract": abstract})

    return records


def attach_query_provenance(abstract_records, concept_from, concept_to):
    """
    Stamp each abstract record with the exact esearch query term that
    surfaced it -- a traceability marker for the write-up and for citation
    review, not an evidence-weighting signal. Under the joint query strategy
    every abstract for a given edge already shares one query (a single AND
    query has no per-record split to carry, unlike separate-then-merge, where
    a record's originating side would actually differ record to record), so
    this is deliberately uniform across all abstracts for the same edge --
    the point is that a reviewer (or the LLM) can see exactly what was
    searched for, right next to each piece of evidence, without re-deriving
    it from the edge's from/to features elsewhere in the output.

    Parameters
    ----------
    abstract_records : list[dict]
        Output of get_abstract_text().
    concept_from, concept_to : str
        The same glossary concept phrases passed to get_relevant_search_results()
        for this edge.

    Returns
    -------
    list[dict]
        abstract_records with a "query_term" key added to each record.
    """
    term = build_search_term(concept_from, concept_to)
    for record in abstract_records:
        record["query_term"] = term
    return abstract_records


# --- INDRA ---

def get_indra_statements(subject=None, obj=None, agents=None, api_key=None,
                          stmt_type=None, max_stmts=25, ev_limit=3):
    """
    Query INDRA for statements linking two agents. Provide EITHER (subject,
    obj) for a directional query -- matches only statement types with
    directional roles (Activation, Inhibition, Phosphorylation,
    IncreaseAmount, DecreaseAmount, etc.) -- OR `agents` (a list of exactly 2
    free-text agent names) for an undirected query, which also catches
    statement types with no subject/object role at all (Complex,
    Translocation, ActiveForm, ...), per the indra_db_service README. Note
    the undirected query is a superset, not a disjoint category: a directional
    statement between the same two agents will also match an undirected
    query, since both agents are still present regardless of role.

    Parameters
    ----------
    subject, obj : str, optional
        Free-text agent names (e.g. "doxorubicin", "neoplasms") or namespace-
        qualified ids (e.g. "TP53@HGNC") for a directional query. INDRA
        resolves plain text via grounding -- works best for literal gene/
        chemical/cytokine names; composite cell-subset/pathway feature names
        are expected to ground to nothing.
    agents : list[str], optional
        Exactly 2 free-text agent names for an undirected query (sent as
        `agent0`/`agent1` -- the REST API's AWS Gateway strips all but one
        `agent` key without the numeric suffix, per the README).
    api_key : str, optional
        INDRA API key. NOT required -- per direct guidance from the INDRA
        team, the from_agents endpoint doesn't need one for our use. Only
        sent as the `api_key` query param if provided.
    stmt_type : str, optional
        Restrict to a specific INDRA statement type, e.g. "Activation".
        Only meaningful for the directional (subject/object) form.
    max_stmts : int
        Max statements to return (INDRA's own cap is 1000).
    ev_limit : int
        Max evidence entries returned per statement.

    Returns
    -------
    list[dict]
        Raw INDRA Statement JSON objects. Empty list if the required
        arguments are missing or the request fails -- never raises.
    """
    if agents:
        if len(agents) != 2 or not agents[0] or not agents[1]:
            return []
        params = {f"agent{i}": agent for i, agent in enumerate(agents)}
    else:
        if not subject or not obj:
            return []
        params = {"subject": subject, "object": obj}
        if stmt_type:
            params["type"] = stmt_type

    params.update({"max_stmts": max_stmts, "ev_limit": ev_limit, "format": "json"})
    if api_key:
        params["api_key"] = api_key

    try:
        response = requests.get(INDRA_STATEMENTS_URL, params=params, timeout=30)
    except requests.exceptions.RequestException:
        return []

    if response.status_code != 200:
        return []

    try:
        payload = response.json()
    except json.JSONDecodeError:
        return []

    return list(payload.get("statements", {}).values())


def _summarize_indra_statements(statements):
    types = collections.Counter(s.get("type") for s in statements)
    pmids = set()
    for statement in statements:
        for evidence in statement.get("evidence", []):
            if evidence.get("pmid"):
                pmids.add(evidence["pmid"])
    return {
        "statement_count": len(statements),
        "types": dict(types),
        "evidence_pmids": sorted(pmids),
    }


def get_directional_indra_support(from_feature, to_feature, api_key=None, max_stmts=25):
    """
    Check INDRA for statements linking from_feature and to_feature, both
    directional (forward = from->to, reverse = to->from) and undirected (any
    statement mentioning both, regardless of role -- catches Complex/
    Translocation/ActiveForm-type statements the directional queries can't
    see; see get_indra_statements() docstring).

    Parameters
    ----------
    from_feature, to_feature : str
        Free-text node names from the bnlearn edge.
    api_key : str, optional
        INDRA API key. Not required -- see get_indra_statements(). Passed
        through only if you have one.
    max_stmts : int
        Passed through to get_indra_statements() for each query.

    Returns
    -------
    dict
        {
          "forward": {"statement_count": int, "types": dict, "evidence_pmids": list},
          "reverse": {"statement_count": int, "types": dict, "evidence_pmids": list},
          "undirected": {"statement_count": int, "types": dict, "evidence_pmids": list},
          "direction_signal": "forward" | "reverse" | "conflicting" | "none",
          "has_any_evidence": bool,
        }
        "direction_signal" is based on forward/reverse only -- an undirected
        hit doesn't claim to support either direction, so it never sets
        direction_signal on its own, but DOES count toward "has_any_evidence".
        Use "has_any_evidence" (not "direction_signal" == "none") to decide
        whether INDRA has ANYTHING on file for this pair.
    """
    forward_statements = get_indra_statements(subject=from_feature, obj=to_feature,
                                               api_key=api_key, max_stmts=max_stmts)
    reverse_statements = get_indra_statements(subject=to_feature, obj=from_feature,
                                               api_key=api_key, max_stmts=max_stmts)
    undirected_statements = get_indra_statements(agents=[from_feature, to_feature],
                                                  api_key=api_key, max_stmts=max_stmts)

    forward_summary = _summarize_indra_statements(forward_statements)
    reverse_summary = _summarize_indra_statements(reverse_statements)
    undirected_summary = _summarize_indra_statements(undirected_statements)

    if forward_summary["statement_count"] and reverse_summary["statement_count"]:
        direction_signal = "conflicting"
    elif forward_summary["statement_count"]:
        direction_signal = "forward"
    elif reverse_summary["statement_count"]:
        direction_signal = "reverse"
    else:
        direction_signal = "none"

    has_any_evidence = bool(
        forward_summary["statement_count"]
        or reverse_summary["statement_count"]
        or undirected_summary["statement_count"]
    )

    return {
        "forward": forward_summary,
        "reverse": reverse_summary,
        "undirected": undirected_summary,
        "direction_signal": direction_signal,
        "has_any_evidence": has_any_evidence,
    }
