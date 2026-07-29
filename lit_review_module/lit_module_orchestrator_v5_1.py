import pandas as pd
import numpy as np

import json
import requests
import os

from lit_module_api_v5_1 import (
    load_feature_glossary,
    get_feature_concept,
    get_relevant_search_results,
    get_abstract_text,
    attach_query_provenance,
    get_directional_indra_support,
)
from lit_module_encoder_v2 import build_corpus, encode_corpus, rank_by_relevance
from lit_module_llm_v5_1 import compare_evidence_to_edge

"""
Orchestrator for the AI-augmented literature review module (v4).

v3 -> v4 change: PubTator3's entity resolver only covers six/seven bioconcept
types and couldn't resolve >50% of this feature set (cell-subset names,
pathway module scores, composite CyTOF readouts, demographic covariates), so
gating the whole pipeline on entity-ID resolution silently skipped the
majority of edges. This version resolves each feature to its plain-English
description from the curated feature glossary (a local lookup, no API call)
and searches PubMed's E-utils (esearch/efetch) with that free text instead --
no entity normalization required.

PubTator itself (both entity resolution and its BioRED relation extraction)
has since been removed from the pipeline entirely -- it only grounds literal
gene/chemical/disease terms, which is a minority of this feature set, and once
INDRA started running unconditionally on every edge (see below), PubTator's
contribution was near-empty signal plus a second external dependency/failure
surface. INDRA is now the pipeline's only structured/typed relationship
source, and it's checked more completely than before: forward, reverse, AND
undirected (statement types like Complex that have no subject/object role and
would otherwise be invisible -- see lit_module_api_v4 docstring).

Pipeline per bnlearn edge (from_feature -> to_feature):
  1. get_feature_concept     : free-text node name -> plain-English glossary
                                description (local lookup, cached in
                                concept_cache; no API call). Required -- a
                                feature missing from the glossary can't be
                                turned into a search query.
  2. get_relevant_search_results : concept descriptions -> list[str] of PMIDs
                                via PubMed esearch (joint from+to AND query --
                                decided over separate-then-merge; see
                                lit_module_api_v4 docstring for rationale). An
                                empty result does NOT skip the row -- PubMed
                                returning nothing is not the same as there
                                being no evidence at all (see below).
  3. get_abstract_text        : PMIDs -> list[dict] with full title/abstract
                                text via PubMed efetch. Skipped (empty list)
                                if step 2 found no PMIDs.
  3b. attach_query_provenance : stamps each abstract record with the exact
                                esearch term that surfaced it -- traceability
                                for the write-up, not an evidence-weighting
                                signal. The edge's output row records the
                                search term regardless of whether it returned
                                anything.
  4. get_directional_indra_support : query INDRA for statements linking the
                                two features -- forward (from->to), reverse
                                (to->from), AND undirected (either role, e.g.
                                Complex formation). Runs on EVERY edge,
                                independent of whether PubMed (steps 2-3)
                                found anything -- INDRA aggregates from
                                curated pathway databases and its own reading
                                pipeline, so it isn't downstream of our
                                esearch call. Also runs without an API key --
                                the INDRA team confirmed one isn't required
                                for this endpoint, so INDRA_API_KEY below is
                                optional, not a gate.

  A row is only marked "inconclusive" and skipped WITHOUT calling the LLM if
  BOTH step 2 (PubMed) found nothing AND step 4's has_any_evidence is False
  (i.e. forward, reverse, AND undirected all came back empty) -- genuinely no
  evidence from either source. This avoids paying for an LLM call with
  nothing behind it, while not discarding real INDRA signal (directional OR
  undirected) just because the PubMed search happened to miss.

  5. build_corpus / encode_corpus : title+abstract -> embeddings (unchanged;
                                returns [] cleanly on an empty abstract list).
  6. rank_by_relevance        : keep the top-k abstracts most similar to the
                                edge (unchanged; returns [] on empty input).
  7. compare_evidence_to_edge : LLM judges whether the literature (PubMed
                                abstracts + INDRA statements) supports the
                                proposed edge AND its direction. Can reach a
                                verdict off INDRA evidence alone, with an
                                empty abstract list, when PubMed found nothing
                                but INDRA did.

Results are written back onto the bnlearn dataframe.
"""

TOP_K_ABSTRACTS = 5

# Per the INDRA team, the from_agents endpoint does NOT require an API key for
# our use (despite the indra_db_service README saying one is required
# elsewhere in the same doc -- confirmed directly with the team instead).
# INDRA_API_KEY is now optional: unset by default, only sent if the team
# issues one later for elevated quota/access.
INDRA_API_KEY = os.environ.get("INDRA_API_KEY")

# Feature glossary (feature_glossary.xlsx) -- override the path via
# FEATURE_GLOSSARY_PATH if it doesn't live next to this script.
FEATURE_GLOSSARY = load_feature_glossary()

# import data from bnlearn output

###### CROSS-VALIDATED INPUTS ###########

# expected columns: "from", "to" (the two nodes of each learned edge)
#full_data = pd.read_csv("./edges/L1b_clinical_curated_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L1a_clinical_full_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L2_olink_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L3a_clinical_full_olink_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L3b_clinical_curated_olink_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L4_wb_pathways_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L5_freq_immunophenotype_literature_edges_cv.csv")
#full_data = pd.read_csv("./edges/L6_pb_signaling_literature_edges_cv.csv")
full_data = pd.read_csv("./edges/L_all_literature_edges_cv.csv")


###### BOOTSTRAP INPUTS ###########

#full_data = pd.read_csv("./boot_edges/L1a_clinical_full_v3_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L1b_clinical_curated_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L2_olink_v2_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L3a_clinical_full_olink_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L3b_clinical_curated_olink_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L4_wb_pathways_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L5_freq_immunophenotype_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L6_pb_signaling_literature_edges.csv")
#full_data = pd.read_csv("./boot_edges/L_all_literature_edges.csv")


in_dag_mask = full_data['in_final_dag'] == True
input_data = full_data[in_dag_mask]


input_data["verdict"] = ""
input_data["confidence_score"] = np.nan
input_data["direction_assessment"] = ""
input_data["rationale"] = ""
# pubmed_references: PMIDs the LLM cited from the PubMed abstracts block only
# (source 1 of the prompt) -- see lit_module_llm_v3's pubmed_citations field.
input_data["pubmed_references"] = ""
input_data["indra_direction_signal"] = ""
input_data["indra_forward_statement_count"] = np.nan
input_data["indra_reverse_statement_count"] = np.nan
input_data["indra_undirected_statement_count"] = np.nan
# indra_*_pmids: written directly from indra_summary's evidence_pmids below,
# NOT from anything the LLM says -- deterministic, so these are never subject
# to whether/how faithfully the model reports what it was shown. Kept
# separate from pubmed_references specifically so it's unambiguous which
# resource (PubMed search vs. INDRA's curated/text-mined statement graph)
# surfaced a given PMID.
input_data["indra_forward_pmids"] = ""
input_data["indra_reverse_pmids"] = ""
input_data["indra_undirected_pmids"] = ""
input_data["search_query"] = ""

# lookup table (limit API calls / repeat lookups -- nodes repeat across edges)
concept_cache = {}


def resolve_feature_concept(feature_name):
    if feature_name not in concept_cache:
        concept_cache[feature_name] = get_feature_concept(feature_name, FEATURE_GLOSSARY)
    return concept_cache[feature_name]


n_total = len(input_data)

print(f"This file has {n_total} edges to process")

for idx, row in input_data.iterrows():

    print(f"working on row {idx+1} of {n_total}")

    from_feature = row["from"]
    to_feature = row["to"]

    # step 1: resolve glossary concept descriptions (local, cached, required)
    concept_from = resolve_feature_concept(from_feature)
    concept_to = resolve_feature_concept(to_feature)

    if not concept_from or not concept_to:
        input_data.at[idx, "verdict"] = "inconclusive"
        input_data.at[idx, "rationale"] = "Feature missing from glossary -- could not build a search query."
        continue

    # step 2: search PubMed (esearch) using the glossary descriptions as free text
    pmids = get_relevant_search_results(concept_from, concept_to)
    input_data.at[idx, "search_query"] = f"({concept_from}) AND ({concept_to})"

    if pmids:
        # step 3: pull full title/abstract text for those PMIDs (efetch)
        abstract_records = get_abstract_text(pmids)

        # step 3b: stamp each abstract with the exact query used, for
        # traceability -- see attach_query_provenance() docstring
        abstract_records = attach_query_provenance(abstract_records, concept_from, concept_to)
    else:
        # No PubMed hits -- NOT the same as "no evidence at all". INDRA is a
        # separate evidence source (curated pathway DBs + its own reading
        # pipeline), not built on top of these PMIDs, so it still gets
        # checked below instead of short-circuiting here.
        abstract_records = []

    # step 4: INDRA -- forward, reverse, and undirected statements. Runs
    # regardless of whether PubMed (steps 2-3) found anything.
    indra_summary = get_directional_indra_support(from_feature, to_feature, INDRA_API_KEY)
    input_data.at[idx, "indra_direction_signal"] = indra_summary["direction_signal"]
    input_data.at[idx, "indra_forward_statement_count"] = indra_summary["forward"]["statement_count"]
    input_data.at[idx, "indra_reverse_statement_count"] = indra_summary["reverse"]["statement_count"]
    input_data.at[idx, "indra_undirected_statement_count"] = indra_summary["undirected"]["statement_count"]

    # PMIDs pulled straight from indra_summary -- deterministic, not routed
    # through the LLM at all, so these are guaranteed complete/accurate and
    # unambiguously INDRA-sourced (see column comments above).
    input_data.at[idx, "indra_forward_pmids"] = ", ".join(indra_summary["forward"]["evidence_pmids"])
    input_data.at[idx, "indra_reverse_pmids"] = ", ".join(indra_summary["reverse"]["evidence_pmids"])
    input_data.at[idx, "indra_undirected_pmids"] = ", ".join(indra_summary["undirected"]["evidence_pmids"])

    if not abstract_records and not indra_summary["has_any_evidence"]:
        # Genuinely nothing from either source -- skip the LLM call rather
        # than paying for a verdict with no evidence behind it.
        input_data.at[idx, "verdict"] = "inconclusive"
        input_data.at[idx, "rationale"] = "No PubMed search results and no INDRA statements for this pair."
        continue

    # steps 5-6: encode and rank by relevance to the proposed edge
    corpus = build_corpus(abstract_records)
    encoded_corpus = encode_corpus(corpus)
    query_text = f"{from_feature} associated with {to_feature}"
    top_abstracts = rank_by_relevance(query_text, encoded_corpus, top_k=TOP_K_ABSTRACTS)

    # step 7: LLM judges the evidence (PubMed abstracts + INDRA statements)
    result = compare_evidence_to_edge(
        from_feature, to_feature, top_abstracts,
        indra_summary=indra_summary,
    )

    input_data.at[idx, "verdict"] = result.get("verdict", "inconclusive")
    input_data.at[idx, "confidence_score"] = result.get("confidence_score", np.nan)
    input_data.at[idx, "direction_assessment"] = result.get("direction_assessment", "no_directional_evidence")
    input_data.at[idx, "rationale"] = result.get("rationale", "")
    # PubMed-only, per the pubmed_citations schema field/prompt instructions
    # in lit_module_llm_v3 -- INDRA's PMIDs live in the indra_*_pmids columns
    # above instead, written directly from indra_summary, not from the model.
    input_data.at[idx, "pubmed_references"] = ", ".join(str(p) for p in result.get("pubmed_citations", []))

############# CROSS-VALIDATED OUTPUTS ###################

#input_data.to_csv("L1b_cv_output_reviewed.csv", index=False)
#input_data.to_csv("L1a_cv_output_reviewed.csv", index=False)
#input_data.to_csv("L2_cv_ouput_reviewed.csv", index=False)
#input_data.to_csv("L3a_cv_ouput_reviewed.csv", index=False)
#input_data.to_csv("L3b_cv_ouput_reviewed.csv", index=False)
#input_data.to_csv("L4_cv_ouput_reviewed.csv", index=False)
#input_data.to_csv("L5_cv_ouput_reviewed.csv", index=False)
#input_data.to_csv("L6_cv_ouput_reviewed.csv", index=False)
input_data.to_csv("Lall_cv_ouput_reviewed.csv", index=False)

############# BOOTSTRAP OUTPUTS ###################

#input_data.to_csv("L1b_boot_output_reviewed.csv", index=False)
#input_data.to_csv("L1a_boot_output_reviewed.csv", index=False)
#input_data.to_csv("L2_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("L3a_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("L3b_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("L4_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("L5_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("L6_boot_ouput_reviewed.csv", index=False)
#input_data.to_csv("Lall_boot_ouput_reviewed.csv", index=False)

concept_cache_df = pd.DataFrame(
    list(concept_cache.items()), columns=["feature_name", "glossary_concept"]
)
concept_cache_df.to_csv("feature_concept_cache.csv", index=False)
