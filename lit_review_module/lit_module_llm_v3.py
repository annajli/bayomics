import pandas as pd
import numpy as np

import os
import json
import requests

"""
LLM comparison module.

Takes the top-ranked, embedding-similarity-filtered abstracts from
encoder_module.rank_by_relevance(), plus two supplemental directionality
signals -- PubTator3's own BioRED relation extractions (summarize_relation_direction)
and INDRA's causal statement graph (get_directional_indra_support) -- and asks an
LLM to judge whether the literature actually supports the bnlearn edge
(from_feature -> to_feature) AND whether the evidence favors that direction
specifically, returning a structured verdict, confidence score, direction
assessment, and citations.

Uses the Anthropic API. Install: pip install anthropic
Set your key: export ANTHROPIC_API_KEY=sk-ant-...
"""

from anthropic import Anthropic

MODEL_NAME = "claude-sonnet-4-5"

_client = None


def _get_client():
    global _client
    if _client is None:
        _client = Anthropic()  # reads ANTHROPIC_API_KEY from env
    return _client


def _format_pubtator_relations(pubtator_relations, from_feature, to_feature):
    if not pubtator_relations:
        return "No PubTator3 relation extractions available for this pair."

    return (
        f'Forward relations ("{from_feature}" -> "{to_feature}"): '
        f"{pubtator_relations.get('forward_count', 0)} found, types: "
        f"{pubtator_relations.get('forward_types') or 'none'}\n"
        f'Reverse relations ("{to_feature}" -> "{from_feature}"): '
        f"{pubtator_relations.get('reverse_count', 0)} found, types: "
        f"{pubtator_relations.get('reverse_types') or 'none'}\n"
        f"PubTator direction signal: {pubtator_relations.get('direction_signal', 'none')}"
    )


def _format_indra_summary(indra_summary, from_feature, to_feature):
    if not indra_summary or indra_summary.get("direction_signal") == "no_api_key":
        return "No INDRA data available (no API key configured)."

    forward = indra_summary.get("forward", {})
    reverse = indra_summary.get("reverse", {})
    return (
        f'Forward statements ("{from_feature}" -> "{to_feature}"): '
        f"{forward.get('statement_count', 0)} found, types: {forward.get('types') or 'none'}\n"
        f'Reverse statements ("{to_feature}" -> "{from_feature}"): '
        f"{reverse.get('statement_count', 0)} found, types: {reverse.get('types') or 'none'}\n"
        f"INDRA direction signal: {indra_summary.get('direction_signal', 'none')}"
    )


def _build_prompt(from_feature, to_feature, ranked_abstracts, pubtator_relations=None, indra_summary=None):
    evidence_blocks = []
    for i, record in enumerate(ranked_abstracts, start=1):
        evidence_blocks.append(
            f"[{i}] PMID {record.get('pmid')} (similarity={record.get('similarity', 0):.3f})\n"
            f"Title: {record.get('title', '')}\n"
            f"Abstract: {record.get('abstract', '')}"
        )
    evidence_text = "\n\n".join(evidence_blocks) if evidence_blocks else "No abstracts retrieved."

    pubtator_text = _format_pubtator_relations(pubtator_relations, from_feature, to_feature)
    indra_text = _format_indra_summary(indra_summary, from_feature, to_feature)

    prompt = f"""You are assisting with an AI-augmented literature review that validates edges
in a Bayesian network learned by bnlearn. The network proposes a DIRECTED relationship:

    "{from_feature}" -> "{to_feature}"

You have three sources of evidence to weigh:

1. PubMed abstracts (ranked by embedding similarity to the proposed relationship):

{evidence_text}

2. PubTator3 BioRED relation extractions (typed, directional, pulled directly
   from the abstracts above):

{pubtator_text}

3. INDRA causal statement graph (aggregates causal/mechanistic statements from
   multiple text-mining sources and curated pathway databases):

{indra_text}

Based on ALL the evidence above, assess (a) whether the literature supports a
relationship between these two entities at all, and (b) whether the evidence
specifically supports the proposed DIRECTION ("{from_feature}" -> "{to_feature}")
rather than the reverse. Weigh sources 2 and 3 most heavily for the direction
question, since they carry explicit directional/causal typing; use source 1 for
general support/context. Respond with ONLY a JSON object, no other text, in this
exact shape:

{{
  "verdict": "supported" | "contradicted" | "inconclusive",
  "confidence_score": <float between 0 and 1>,
  "direction_assessment": "supports_proposed_direction" | "supports_reverse_direction" | "conflicting" | "no_directional_evidence",
  "rationale": "<1-3 sentence explanation>",
  "citations": [<list of PMIDs actually used as evidence>]
}}"""
    return prompt


def compare_evidence_to_edge(from_feature, to_feature, ranked_abstracts, pubtator_relations=None, indra_summary=None):
    """
    Ask the LLM to judge whether the retrieved literature supports a bnlearn edge,
    and specifically whether it supports the proposed direction.

    Parameters
    ----------
    from_feature, to_feature : str
        The two node names from the bnlearn edge being evaluated.
    ranked_abstracts : list[dict]
        Output of encoder_module.rank_by_relevance() -- records with "pmid",
        "title", "abstract", "similarity".
    pubtator_relations : dict, optional
        Output of pubtator_api_module.summarize_relation_direction(). If omitted,
        the prompt notes no PubTator relation data was supplied.
    indra_summary : dict, optional
        Output of pubtator_api_module.get_directional_indra_support(). If omitted,
        the prompt notes no INDRA data was supplied.

    Returns
    -------
    dict
        {"verdict": str, "confidence_score": float, "direction_assessment": str,
         "rationale": str, "citations": list}
        Falls back to a default "inconclusive" dict if parsing fails.
    """
    client = _get_client()
    prompt = _build_prompt(from_feature, to_feature, ranked_abstracts, pubtator_relations, indra_summary)

    response = client.messages.create(
        model=MODEL_NAME,
        max_tokens=500,
        messages=[{"role": "user", "content": prompt}],
    )

    raw_text = response.content[0].text.strip()

    try:
        result = json.loads(raw_text)
    except json.JSONDecodeError:
        result = {
            "verdict": "inconclusive",
            "confidence_score": 0.0,
            "direction_assessment": "no_directional_evidence",
            "rationale": f"Could not parse LLM response: {raw_text[:200]}",
            "citations": [],
        }

    return result