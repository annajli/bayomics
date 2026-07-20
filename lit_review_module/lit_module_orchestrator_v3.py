import pandas as pd
import numpy as np

import json
import requests
import os

from pubtator_api_module_v3 import (
    get_entity_ID,
    get_relevant_search_results,
    get_abstract_text,
    summarize_relation_direction,
    get_directional_indra_support,
)
from lit_module_encoder_v2 import build_corpus, encode_corpus, rank_by_relevance
from lit_module_llm_v3 import compare_evidence_to_edge

"""
Orchestrator for the AI-augmented literature review module.

Pipeline per bnlearn edge (from_feature -> to_feature):
  1. get_entity_ID          : free-text node name -> PubTator3 accession ID
                               (cached in entity_id_dict so repeated nodes only
                               cost one API call each)
  2. get_relevant_search_results : accession IDs -> list[dict] of candidate papers
                               (title/snippet/metadata only, no full abstract)
  3. get_abstract_text       : PMIDs from step 2 -> list[dict] with full abstract
                               text AND raw BioRED relations (the /search/
                               endpoint alone doesn't give you enough text to
                               encode meaningfully, nor any relation typing)
  4. summarize_relation_direction : classify PubTator's own relation extractions
                               as forward/reverse/conflicting/none relative to
                               the proposed edge direction -- free signal, no
                               extra API call
  5. get_directional_indra_support : query INDRA for causal statements in both
                               the proposed direction and the reverse, to
                               supplement PubTator (which has no dedicated
                               causal-direction relation type) with an explicit
                               subject->object causal graph check
  6. build_corpus / encode_corpus : title+abstract -> embeddings
  7. rank_by_relevance       : keep the top-k abstracts most similar to the edge
  8. compare_evidence_to_edge : LLM judges whether the literature (abstracts +
                               PubTator relations + INDRA statements) supports
                               the proposed edge AND its direction, returns a
                               verdict + confidence_score + direction_assessment
                               + citations

Results are written back onto the bnlearn dataframe.
"""

TOP_K_ABSTRACTS = 5

# INDRA requires an API key (request one from the INDRA team). If unset, INDRA
# lookups are skipped and direction_assessment falls back to PubTator + LLM only.
INDRA_API_KEY = os.environ.get("INDRA_API_KEY")
#INDRA_API_KEY = os.environ.get("indra_api_key")

# import data from bnlearn output
# expected columns: "from", "to" (the two nodes of each learned edge)
input_data = pd.read_csv("bnlearn_output.csv")
input_data["verdict"] = ""
input_data["confidence_score"] = np.nan
input_data["direction_assessment"] = ""
input_data["rationale"] = ""
input_data["references"] = ""
input_data["pubtator_direction_signal"] = ""
input_data["pubtator_forward_relation_types"] = ""
input_data["pubtator_reverse_relation_types"] = ""
input_data["indra_direction_signal"] = ""
input_data["indra_forward_statement_count"] = np.nan
input_data["indra_reverse_statement_count"] = np.nan

# make a lookup table for entity IDs (limit API calls -- nodes repeat across edges)
entity_id_dict = {}


def resolve_entity_id(feature_name):
    if feature_name not in entity_id_dict:
        entity_id_dict[feature_name] = get_entity_ID(feature_name)
    return entity_id_dict[feature_name]


for idx, row in input_data.iterrows():
    from_feature = row["from"]
    to_feature = row["to"]

    # step 1: resolve entity IDs (cached)
    from_id = resolve_entity_id(from_feature)
    to_id = resolve_entity_id(to_feature)

    if not from_id or not to_id:
        input_data.at[idx, "verdict"] = "inconclusive"
        input_data.at[idx, "rationale"] = "Could not resolve one or both entity IDs."
        continue

    # step 2: search for relevant papers (list[dict], title/snippet only)
    search_results = get_relevant_search_results(from_id, to_id)

    if not search_results:
        input_data.at[idx, "verdict"] = "inconclusive"
        input_data.at[idx, "rationale"] = "No PubTator3 search results for this pair."
        continue

    # step 3: pull full abstract text (+ raw BioRED relations) for those papers
    pmids = [r["pmid"] for r in search_results if r.get("pmid")]
    abstract_records = get_abstract_text(pmids)

    # step 4: free directional signal from PubTator's own relation extractions
    pubtator_relations = summarize_relation_direction(abstract_records, from_id, to_id)
    input_data.at[idx, "pubtator_direction_signal"] = pubtator_relations["direction_signal"]
    input_data.at[idx, "pubtator_forward_relation_types"] = ", ".join(pubtator_relations["forward_types"])
    input_data.at[idx, "pubtator_reverse_relation_types"] = ", ".join(pubtator_relations["reverse_types"])

    # step 5: supplemental directional signal from INDRA's causal statement graph
    indra_summary = get_directional_indra_support(from_feature, to_feature, INDRA_API_KEY)
    input_data.at[idx, "indra_direction_signal"] = indra_summary["direction_signal"]
    input_data.at[idx, "indra_forward_statement_count"] = indra_summary["forward"]["statement_count"]
    input_data.at[idx, "indra_reverse_statement_count"] = indra_summary["reverse"]["statement_count"]

    # steps 6-7: encode and rank by relevance to the proposed edge
    corpus = build_corpus(abstract_records)
    encoded_corpus = encode_corpus(corpus)
    query_text = f"{from_feature} associated with {to_feature}"
    top_abstracts = rank_by_relevance(query_text, encoded_corpus, top_k=TOP_K_ABSTRACTS)

    # step 8: LLM judges the evidence (abstracts + PubTator relations + INDRA statements)
    result = compare_evidence_to_edge(
        from_feature, to_feature, top_abstracts,
        pubtator_relations=pubtator_relations,
        indra_summary=indra_summary,
    )

    input_data.at[idx, "verdict"] = result.get("verdict", "inconclusive")
    input_data.at[idx, "confidence_score"] = result.get("confidence_score", np.nan)
    input_data.at[idx, "direction_assessment"] = result.get("direction_assessment", "no_directional_evidence")
    input_data.at[idx, "rationale"] = result.get("rationale", "")
    input_data.at[idx, "references"] = ", ".join(str(p) for p in result.get("citations", []))

input_data.to_csv("bnlearn_output_reviewed.csv", index=False)