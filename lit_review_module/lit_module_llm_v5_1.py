import pandas as pd
import numpy as np

import os
import requests

"""
LLM comparison module.

Takes the top-ranked, embedding-similarity-filtered abstracts from
encoder_module.rank_by_relevance(), plus INDRA's causal statement graph
(get_directional_indra_support) -- and asks an LLM to judge whether the
literature actually supports the bnlearn edge (from_feature -> to_feature)
AND whether the evidence favors that direction specifically, returning a
structured verdict, confidence score, direction assessment, and PubMed
citations (INDRA's PMIDs are tracked separately -- see pubmed_citations note
on compare_evidence_to_edge below, and lit_module_orchestrator_v4 for the
INDRA-side PMID columns).

PubTator's BioRED relation extractions were dropped as a second structured
signal (see lit_module_api_v4 module docstring for why): PubTator only
grounds literal gene/chemical/disease terms, so on this feature set it
mostly returned empty relations and added noise/failure surface rather than
real evidence, now that INDRA runs unconditionally on every edge (including
an undirected query, so non-directional relationship types aren't missed
either). INDRA is now the pipeline's only structured/typed evidence source,
alongside the free-text PubMed abstracts.

Verdict is returned via a FORCED TOOL CALL (record_verdict, defined below)
rather than free-text JSON parsing. The original approach asked the model to
"respond with ONLY a JSON object" and ran json.loads() on the raw text, which
broke in production two ways: Claude would sometimes wrap the object in a
markdown code fence (```json ... ```) despite the instruction not to, and long
rationales occasionally pushed the response past the old max_tokens=500 cap,
truncating it mid-object. Both produced the same symptom -- a
json.JSONDecodeError caught by the except block, surfacing as "Could not parse
LLM response: ..." in the rationale column. Forcing a tool call sidesteps
both: Anthropic enforces the input_schema server-side, so there's no
markdown/prose wrapping possible, and no text parsing step at all --
response.content contains a tool_use block whose .input is already a
schema-conformant dict.

Uses the Anthropic API. Install: pip install anthropic
Set your key: export ANTHROPIC_API_KEY=sk-ant-...
"""

from anthropic import Anthropic
from dotenv import load_dotenv



MODEL_NAME = "claude-sonnet-4-5"

_client = None

load_dotenv()

if not os.getenv("ANTHROPIC_API_KEY"):
    raise RuntimeError("ANTHROPIC_API_KEY not set — check your .env file")

def _get_client():
    global _client
    if _client is None:
        _client = Anthropic()  # reads ANTHROPIC_API_KEY from env
    return _client


TOOL_NAME = "record_verdict"

VERDICT_TOOL = {
    "name": TOOL_NAME,
    "description": (
        "Record the structured assessment of whether the retrieved literature "
        "supports the proposed bnlearn edge and, specifically, whether the "
        "evidence favors the proposed direction over the reverse."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "verdict": {
                "type": "string",
                "enum": ["supported", "contradicted", "inconclusive"],
                "description": "Whether the literature supports a relationship between the two entities at all.",
            },
            "confidence_score": {
                "type": "number",
                "minimum": 0,
                "maximum": 1,
                "description": "Confidence in the verdict, from 0 (none) to 1 (fully confident).",
            },
            "direction_assessment": {
                "type": "string",
                "enum": [
                    "supports_proposed_direction",
                    "supports_reverse_direction",
                    "conflicting",
                    "no_directional_evidence",
                ],
                "description": "Whether the evidence specifically supports the proposed from->to direction, the reverse, both (conflicting), or neither.",
            },
            "rationale": {
                "type": "string",
                "description": "1-3 sentence explanation for the verdict and direction assessment.",
            },
            "pubmed_citations": {
                "type": "array",
                "items": {"type": "string"},
                "description": "PMIDs from the PubMed abstracts block (source 1) ONLY that were actually used as evidence for this verdict. Do NOT include PMIDs from the INDRA statements block (source 2) here, even though some are listed there too -- those are tracked separately by the pipeline, not by this field.",
            },
        },
        "required": ["verdict", "confidence_score", "direction_assessment", "rationale", "pubmed_citations"],
    },
}

# Cap how many PMIDs get listed per INDRA bucket (forward/reverse/undirected)
# in the prompt -- get_indra_statements can return up to max_stmts statements
# with up to ev_limit evidence entries each, so evidence_pmids can run long;
# uncapped, that risks bloating every prompt for a handful of well-studied
# gene pairs. Statement/type counts above are never truncated, only this list.
INDRA_PMID_DISPLAY_LIMIT = 15


def _format_pmid_list(pmids, limit=INDRA_PMID_DISPLAY_LIMIT):
    if not pmids:
        return "none"
    shown = ", ".join(pmids[:limit])
    if len(pmids) > limit:
        shown += f", ... ({len(pmids) - limit} more not shown)"
    return shown


def _format_indra_summary(indra_summary, from_feature, to_feature):
    if not indra_summary:
        return "No INDRA data available for this pair."

    forward = indra_summary.get("forward", {})
    reverse = indra_summary.get("reverse", {})
    undirected = indra_summary.get("undirected", {})
    return (
        f'Forward statements ("{from_feature}" -> "{to_feature}"): '
        f"{forward.get('statement_count', 0)} found, types: {forward.get('types') or 'none'}, "
        f"PMIDs: {_format_pmid_list(forward.get('evidence_pmids'))}\n"
        f'Reverse statements ("{to_feature}" -> "{from_feature}"): '
        f"{reverse.get('statement_count', 0)} found, types: {reverse.get('types') or 'none'}, "
        f"PMIDs: {_format_pmid_list(reverse.get('evidence_pmids'))}\n"
        f"Undirected statements (either role -- e.g. complex formation or other "
        f"non-directional relationship types not captured above; does NOT by "
        f"itself support either direction): "
        f"{undirected.get('statement_count', 0)} found, types: {undirected.get('types') or 'none'}, "
        f"PMIDs: {_format_pmid_list(undirected.get('evidence_pmids'))}\n"
        f"INDRA direction signal: {indra_summary.get('direction_signal', 'none')}"
    )


def _build_prompt(from_feature, to_feature, ranked_abstracts, indra_summary=None):
    evidence_blocks = []
    for i, record in enumerate(ranked_abstracts, start=1):
        # query_term (added in lit_module_api_v4.attach_query_provenance) is a
        # traceability marker, not a weighting signal -- shown here so the LLM
        # (and a human reviewer reading the rationale later) can see exactly
        # what was searched for, without needing to cross-reference the edge's
        # from/to features elsewhere in the output.
        query_line = f"Retrieved via query: {record['query_term']}\n" if record.get("query_term") else ""
        evidence_blocks.append(
            f"[{i}] PMID {record.get('pmid')} (similarity={record.get('similarity', 0):.3f})\n"
            f"{query_line}"
            f"Title: {record.get('title', '')}\n"
            f"Abstract: {record.get('abstract', '')}"
        )
    evidence_text = "\n\n".join(evidence_blocks) if evidence_blocks else "No abstracts retrieved."

    indra_text = _format_indra_summary(indra_summary, from_feature, to_feature)

    prompt = f"""You are assisting with an AI-augmented literature review that validates edges
in a Bayesian network learned by bnlearn. The network proposes a DIRECTED relationship:

    "{from_feature}" -> "{to_feature}"

You have two sources of evidence to weigh:

1. PubMed abstracts (ranked by embedding similarity to the proposed relationship):

{evidence_text}

2. INDRA causal statement graph (aggregates causal/mechanistic statements from
   multiple text-mining sources and curated pathway databases; includes both
   directional statements, e.g. Activation/Inhibition, and undirected
   statements, e.g. Complex formation, that don't imply a direction):

{indra_text}

Based on ALL the evidence above, assess (a) whether the literature supports a
relationship between these two entities at all, and (b) whether the evidence
specifically supports the proposed DIRECTION ("{from_feature}" -> "{to_feature}")
rather than the reverse. Weigh source 2's forward/reverse counts most heavily
for the direction question, since they carry explicit directional/causal
typing; an INDRA "undirected" count supports that a relationship exists but
should NOT be treated as evidence for either direction. Use source 1 for
general support/context and for direction only when the abstract text itself
states a direction. For pubmed_citations, list ONLY PMIDs from source 1 (the
PubMed abstracts) -- the pipeline tracks INDRA's PMIDs separately, so do not
put source 2's PMIDs in this field even though your reasoning may rely on
them. Call the {TOOL_NAME} tool with your assessment."""
    return prompt


def compare_evidence_to_edge(from_feature, to_feature, ranked_abstracts, indra_summary=None):
    """
    Ask the LLM to judge whether the retrieved literature supports a bnlearn edge,
    and specifically whether it supports the proposed direction.

    Parameters
    ----------
    from_feature, to_feature : str
        The two node names from the bnlearn edge being evaluated.
    ranked_abstracts : list[dict]
        Output of encoder_module.rank_by_relevance() -- records with "pmid",
        "title", "abstract", "similarity", and "query_term" -- the exact
        esearch query that surfaced the record, shown in the prompt for
        traceability, not used as a weighting signal.
    indra_summary : dict, optional
        Output of lit_module_api_v4.get_directional_indra_support() -- has
        "forward", "reverse", "undirected", "direction_signal", and
        "has_any_evidence" keys. Each of forward/reverse/undirected carries an
        "evidence_pmids" list, which IS shown to the model (capped at
        INDRA_PMID_DISPLAY_LIMIT per bucket) for the LLM's own reasoning, but
        is intentionally excluded from "pubmed_citations" below -- INDRA's
        PMIDs are tracked deterministically by the orchestrator directly from
        indra_summary (not from anything the model says), precisely so the
        two sources never get commingled in one ambiguous citation list. If
        indra_summary is omitted, the prompt notes no INDRA data was
        supplied.

    Returns
    -------
    dict
        {"verdict": str, "confidence_score": float, "direction_assessment": str,
         "rationale": str, "pubmed_citations": list}
        "pubmed_citations" is PubMed-abstract PMIDs ONLY, per the tool schema
        and prompt instructions -- see lit_module_orchestrator_v4 for where
        INDRA's PMIDs are separately written to their own output columns.
        Falls back to a default "inconclusive" dict only if the model
        somehow declines the forced tool call (e.g. a content refusal) --
        no longer falls back on ordinary JSON-formatting variation, since
        there's no free-text JSON to format correctly in the first place.
    """
    client = _get_client()
    prompt = _build_prompt(from_feature, to_feature, ranked_abstracts, indra_summary)

    response = client.messages.create(
        model=MODEL_NAME,
        max_tokens=800,
        tools=[VERDICT_TOOL],
        tool_choice={"type": "tool", "name": TOOL_NAME},
        messages=[{"role": "user", "content": prompt}],
    )

    tool_use_block = next(
        (block for block in response.content if getattr(block, "type", None) == "tool_use"),
        None,
    )

    if tool_use_block is None:
        # Shouldn't normally happen with tool_choice forcing the call, but
        # guard against a refusal/empty response rather than crashing the run.
        fallback_text = "".join(
            getattr(block, "text", "") for block in response.content if getattr(block, "type", None) == "text"
        )
        return {
            "verdict": "inconclusive",
            "confidence_score": 0.0,
            "direction_assessment": "no_directional_evidence",
            "rationale": f"Model did not return a tool call (stop_reason={response.stop_reason}): {fallback_text[:200]}",
            "pubmed_citations": [],
        }

    payload = tool_use_block.input
    return {
        "verdict": payload.get("verdict", "inconclusive"),
        "confidence_score": payload.get("confidence_score", 0.0),
        "direction_assessment": payload.get("direction_assessment", "no_directional_evidence"),
        "rationale": payload.get("rationale", ""),
        "pubmed_citations": payload.get("pubmed_citations", []),
    }
