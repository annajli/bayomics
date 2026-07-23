import pandas as pd
import streamlit as st

from data.loaders import load_edges
from viz.graph import filter_edges

ROOTS = {
    "age_group": ("Age", "→ Olink proteins + cell frequencies",
                  "Older adults show elevated GDF15 and declining naive T-cell fractions."),
    "sex": ("Biological sex", "→ Olink proteins + clinical labs",
            "Sex drives LEP, IGFBP7, and creatinine — no cell-frequency or signaling effect."),
    "cmv": ("CMV serostatus", "→ Cell composition only",
            "CMV reshapes cell frequencies (cytotoxic memory CD4, adaptive NK, γδ T "
            "effectors) but does not drive proteins, labs, or per-cell signaling."),
}


def _edges_from_root(root: str) -> pd.DataFrame:
    """All final-DAG edges originating at a root, across every network, deduped."""
    frames = []
    from config.network_registry import LEVEL_ORDER
    for lvl in LEVEL_ORDER:
        try:
            e = filter_edges(load_edges(lvl), threshold=0.0, dag_only=True)
        except FileNotFoundError:
            continue
        e = e[e["from"] == root][["to", "strength"]].copy()
        e["level"] = lvl
        frames.append(e)
    if not frames:
        return pd.DataFrame(columns=["to", "strength", "level"])
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values("strength", ascending=False)


def render():
    st.title("🧭 Root → Modality Findings")
    st.markdown(
        "Three discrete roots — **age, sex, CMV** — anchor every network. Their downstream "
        "effects are **modality-specific and robust** across the ladder, including the "
        "352-node L_all model. Edges below are the direct root→feature links present in the "
        "averaged DAGs, pulled live from the edge lists."
    )

    for root, (name, arrow, blurb) in ROOTS.items():
        st.subheader(f"{name}  {arrow}")
        st.caption(blurb)
        df = _edges_from_root(root)
        if len(df) == 0:
            st.info("No direct edges from this root in any final DAG.")
        else:
            st.dataframe(
                df.rename(columns={"to": "Target feature", "strength": "Strength",
                                   "level": "Network"}),
                hide_index=True, use_container_width=True,
                height=min(300, 40 + 35 * len(df)),
            )
        st.divider()
