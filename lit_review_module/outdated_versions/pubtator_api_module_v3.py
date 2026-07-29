import pandas as pd
import numpy as np

import os
import json
import time
import collections
import requests

"""
PubTator3 + INDRA API wrapper functions.

--- PubTator3 ---

1) Find Entity ID (autocomplete):
   https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/?query=QUERY&concept=BIOCONCEPT(OPTIONAL)&limit=NUMBER(OPTIONAL)

   Returns a JSON list of candidate entities, e.g.:
   [{"_id": "@CHEMICAL_Doxorubicin", "biotype": "chemical", "db_id": "D004317",
     "db": "ncbi_mesh", "name": "Doxorubicin", "match": "Matched on name <m>Doxorubicin</m>"}]

   The "_id" field (e.g. "@CHEMICAL_Doxorubicin") is the accession string PubTator3
   expects when you reference this entity in a search/relation query later.

2) Export relevant search results:
   https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text=QUERY

   Returns a JSON object (dict) shaped like:
   {"results": [ {...paper 1...}, {...paper 2...}, ... ],
    "facets": {...}, "page_size": 10, "current": 1, "count": N, "total_pages": M}

   Each item under "results" is a dict with keys like pmid, title, journal, authors,
   score, text_hl (a short highlighted snippet -- NOT the full abstract), citations.
   There is no full abstract text in this response.

3) Fetch full abstract text + relations (BioC-JSON export):
   https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids=PMID1,PMID2,...

   Returns {"PubTator3": [ {...article with "passages" and "relations"...}, ... ]}.
   Each article has a "passages" list; the passage with infons.type == "abstract"
   holds the full abstract text in passage["text"]. This is the actual text you
   want to hand to the encoder module.

   Each article ALSO has a top-level "relations" list -- PubTator3's own BioRED
   relation extractor output. Each relation is typed and directional:
   {"infons": {"type": "Negative_Correlation", "score": "0.9988",
               "role1": {"name": "Doxorubicin", "accession": "@CHEMICAL_Doxorubicin"},
               "role2": {"name": "Neoplasms", "accession": "@DISEASE_Neoplasms"}}}
   role1 -> role2 gives a directional signal for free, at no extra API cost --
   see summarize_relation_direction() below.

--- INDRA ---

INDRA (Integrated Network and Dynamical Reasoning Assembler) aggregates causal /
mechanistic statements from multiple text-mining readers and curated pathway
databases into a single typed format with explicit subject/object roles, e.g.
Activation, Inhibition, IncreaseAmount, DecreaseAmount. Querying subject=X&object=Y
directly tests whether the literature/knowledge-graph supports THAT direction --
a more direct directionality check than PubTator's symmetric "associate" relation.

    GET https://db.indra.bio/statements/from_agents?subject=X&object=Y&api_key=KEY

Requires an API key (request one from the INDRA team). See
get_indra_statements() and get_directional_indra_support() below.
"""

AUTOCOMPLETE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/"
SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
EXPORT_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
INDRA_STATEMENTS_URL = "https://db.indra.bio/statements/from_agents"


def get_entity_ID(input_feature, concept=None, limit=1):
    """
    Look up PubTator3's known entity ID (accession string) for a free-text feature
    name coming out of the bnlearn network (e.g. "doxorubicin").

    Parameters
    ----------
    input_feature : str
        Free-text entity name, e.g. "doxorubicin".
    concept : str, optional
        Restrict to a bioconcept type (e.g. "chemical", "disease", "gene").
    limit : int
        Number of candidates to request from the API (we only use the first).

    Returns
    -------
    str or None
        The accession string, e.g. "@CHEMICAL_Doxorubicin", or None if no match
        was found / the request failed.
    """
    params = {"query": input_feature, "limit": limit}
    if concept:
        params["concept"] = concept

    response = requests.get(AUTOCOMPLETE_URL, params=params, timeout=15)

    if response.status_code != 200:
        return None

    candidates = response.json()  # this endpoint returns a JSON list
    if not candidates:
        return None

    return candidates[0]["_id"]


def get_relevant_search_results(id_from, id_to, relationship="associate", max_pages=1):
    """
    Query PubTator3 for publications relating two known entity IDs.

    Parameters
    ----------
    id_from, id_to : str
        Accession strings returned by get_entity_ID, e.g. "@CHEMICAL_Doxorubicin".
    relationship : str
        Relationship qualifier PubTator3 accepts (e.g. "associate", "treat", "ANY").
    max_pages : int
        How many result pages to pull (10 results/page). 1 is usually enough for
        a lit-review pass; raise this if you want a deeper search per edge.

    Returns
    -------
    list[dict]
        A list of paper-record dictionaries (the "results" array from the API
        response). Each dict includes pmid, title, journal, authors, score, and
        a highlighted snippet (text_hl) -- but NOT the full abstract. Returns an
        empty list if either ID is missing or the request fails.
    """
    if not id_from or not id_to:
        return []

    text_param = f"relations:{relationship}|{id_from}|{id_to}"
    all_results = []

    for page in range(1, max_pages + 1):
        response = requests.get(
            SEARCH_URL, params={"text": text_param, "page": page}, timeout=15
        )
        if response.status_code != 200:
            break

        payload = response.json()  # this endpoint returns a JSON object (dict)
        page_results = payload.get("results", [])
        all_results.extend(page_results)

        if page >= payload.get("total_pages", 1):
            break

    return all_results


def get_abstract_text(pmids, batch_size=100):
    """
    Fetch full abstract text for a list of PMIDs via the BioC-JSON export endpoint.
    This is the piece missing from the original two-call flow: the /search/ endpoint
    only returns a short snippet, not the full abstract needed for encoding.

    Parameters
    ----------
    pmids : list[int|str] or int|str
        One PMID or a list of PMIDs (e.g. pulled from get_relevant_search_results
        output via [r["pmid"] for r in results]).
    batch_size : int
        Number of PMIDs to request per call (keeps URLs/requests reasonably sized).

    Returns
    -------
    list[dict]
        One dict per publication: {"pmid": ..., "title": ..., "abstract": ...,
        "relations": [...]}. "abstract" is "" if no abstract passage was found.
        "relations" is the article's raw BioRED relation list -- each item is
        {"type": str, "score": float or None,
         "role1_name": str, "role1_id": str,
         "role2_name": str, "role2_id": str}. See summarize_relation_direction().
    """
    if pmids is None:
        return []
    if not isinstance(pmids, (list, tuple, set)):
        pmids = [pmids]
    pmids = [str(p) for p in pmids]

    records = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        response = requests.get(
            EXPORT_URL, params={"pmids": ",".join(batch)}, timeout=30
        )
        if response.status_code != 200:
            continue

        payload = response.json()
        for article in payload.get("PubTator3", []):
            title = ""
            abstract = ""
            for passage in article.get("passages", []):
                passage_type = passage.get("infons", {}).get("type")
                if passage_type == "title":
                    title = passage.get("text", "")
                elif passage_type == "abstract":
                    abstract = passage.get("text", "")

            relations = []
            for relation in article.get("relations", []):
                infons = relation.get("infons", {})
                role1 = infons.get("role1") or {}
                role2 = infons.get("role2") or {}
                score = infons.get("score")
                try:
                    score = float(score) if score is not None else None
                except (TypeError, ValueError):
                    score = None
                relations.append(
                    {
                        "type": infons.get("type"),
                        "score": score,
                        "role1_name": role1.get("name"),
                        "role1_id": role1.get("accession"),
                        "role2_name": role2.get("name"),
                        "role2_id": role2.get("accession"),
                    }
                )

            records.append(
                {
                    "pmid": article.get("pmid"),
                    "title": title,
                    "abstract": abstract,
                    "relations": relations,
                }
            )

    return records


def summarize_relation_direction(abstract_records, from_id, to_id):
    """
    Classify the directionality of PubTator3's own BioRED relation extractions
    (already returned by get_abstract_text) relative to a proposed bnlearn edge
    from_id -> to_id. Free directional signal -- no extra API call.

    Parameters
    ----------
    abstract_records : list[dict]
        Output of get_abstract_text() (each record has a "relations" list).
    from_id, to_id : str
        Accession strings for the proposed edge, e.g. "@CHEMICAL_Doxorubicin",
        "@DISEASE_Neoplasms".

    Returns
    -------
    dict
        {
          "forward_count": int,   # relations where role1==from_id, role2==to_id
          "reverse_count": int,   # relations where role1==to_id, role2==from_id
          "forward_types": list[str],  # distinct relation types seen forward
          "reverse_types": list[str],  # distinct relation types seen reverse
          "direction_signal": "forward" | "reverse" | "conflicting" | "none",
        }
    """
    forward_types = collections.Counter()
    reverse_types = collections.Counter()

    for record in abstract_records:
        for relation in record.get("relations", []):
            rel_type = relation.get("type")
            if relation.get("role1_id") == from_id and relation.get("role2_id") == to_id:
                forward_types[rel_type] += 1
            elif relation.get("role1_id") == to_id and relation.get("role2_id") == from_id:
                reverse_types[rel_type] += 1

    forward_count = sum(forward_types.values())
    reverse_count = sum(reverse_types.values())

    if forward_count and reverse_count:
        direction_signal = "conflicting"
    elif forward_count:
        direction_signal = "forward"
    elif reverse_count:
        direction_signal = "reverse"
    else:
        direction_signal = "none"

    return {
        "forward_count": forward_count,
        "reverse_count": reverse_count,
        "forward_types": sorted(forward_types.keys()),
        "reverse_types": sorted(reverse_types.keys()),
        "direction_signal": direction_signal,
    }


def get_indra_statements(subject, obj, api_key, stmt_type=None, max_stmts=25, ev_limit=3):
    """
    Query INDRA for causal statements with a specific subject -> object direction.
    Only statements matching THIS direction are returned (INDRA does not return
    the reverse direction from a single call) -- see get_directional_indra_support()
    for a combined forward + reverse check.

    Parameters
    ----------
    subject, obj : str
        Free-text agent names (e.g. "doxorubicin", "neoplasms") or namespace-
        qualified ids (e.g. "TP53@HGNC"). INDRA resolves plain text via grounding.
    api_key : str
        INDRA API key (request one from the INDRA team).
    stmt_type : str, optional
        Restrict to a specific INDRA statement type, e.g. "Activation", "Inhibition".
    max_stmts : int
        Max statements to return (INDRA's own cap is 1000).
    ev_limit : int
        Max evidence entries returned per statement.

    Returns
    -------
    list[dict]
        Raw INDRA Statement JSON objects (each has "type" and an "evidence" list
        with "pmid"/"text"). Returns an empty list if the API key is missing or
        the request fails.
    """
    if not api_key or not subject or not obj:
        return []

    params = {
        "subject": subject,
        "object": obj,
        "api_key": api_key,
        "max_stmts": max_stmts,
        "ev_limit": ev_limit,
        "format": "json",
    }
    if stmt_type:
        params["type"] = stmt_type

    response = requests.get(INDRA_STATEMENTS_URL, params=params, timeout=30)
    if response.status_code != 200:
        return []

    payload = response.json()
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


def get_directional_indra_support(from_feature, to_feature, api_key, max_stmts=25):
    """
    Check INDRA for causal statements in both the proposed bnlearn direction
    (from_feature -> to_feature) and the reverse, so you can see whether the
    literature/knowledge-graph evidence actually favors the direction bnlearn
    proposed, the opposite direction, both, or neither.

    Parameters
    ----------
    from_feature, to_feature : str
        Free-text node names from the bnlearn edge.
    api_key : str
        INDRA API key. If falsy, returns a "no_api_key" result without calling out.
    max_stmts : int
        Passed through to get_indra_statements() for each direction.

    Returns
    -------
    dict
        {
          "forward": {"statement_count": int, "types": dict, "evidence_pmids": list},
          "reverse": {"statement_count": int, "types": dict, "evidence_pmids": list},
          "direction_signal": "forward" | "reverse" | "conflicting" | "none" | "no_api_key",
        }
    """
    if not api_key:
        empty = {"statement_count": 0, "types": {}, "evidence_pmids": []}
        return {"forward": empty, "reverse": empty, "direction_signal": "no_api_key"}

    forward_statements = get_indra_statements(from_feature, to_feature, api_key, max_stmts=max_stmts)
    reverse_statements = get_indra_statements(to_feature, from_feature, api_key, max_stmts=max_stmts)

    forward_summary = _summarize_indra_statements(forward_statements)
    reverse_summary = _summarize_indra_statements(reverse_statements)

    if forward_summary["statement_count"] and reverse_summary["statement_count"]:
        direction_signal = "conflicting"
    elif forward_summary["statement_count"]:
        direction_signal = "forward"
    elif reverse_summary["statement_count"]:
        direction_signal = "reverse"
    else:
        direction_signal = "none"

    return {
        "forward": forward_summary,
        "reverse": reverse_summary,
        "direction_signal": direction_signal,
    }