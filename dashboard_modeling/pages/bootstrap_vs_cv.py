import pandas as pd
import streamlit as st

from config.network_registry import NETWORKS, LEVEL_ORDER
from data.loaders import load_edges
from viz.graph import filter_edges


def _edge_pairs(df: pd.DataFrame) -> set:
    return {frozenset((r["from"], r["to"])) for _, r in df.iterrows()}


def _pretty(pairs: set) -> pd.DataFrame:
    rows = [sorted(list(p)) for p in pairs]
    return pd.DataFrame(rows, columns=["Node A", "Node B"]).sort_values("Node A")


def render():
    st.title("⚖️ Bootstrap vs Cross-Validation")
    st.markdown(
        "The **bootstrap** pipeline (`bn_learning/`) measures edge stability under "
        "resampling; the **cross-validation** pipeline (`bn_cv/`) measures predictive "
        "generalization. Edges present in *both* final DAGs are the most trustworthy."
    )

    level = st.selectbox("Network level", LEVEL_ORDER,
                         format_func=lambda lvl: NETWORKS[lvl]["label"])

    try:
        boot = filter_edges(load_edges(level, "bootstrap"), 0.0, dag_only=True)
        cv = filter_edges(load_edges(level, "cv"), 0.0, dag_only=True)
    except FileNotFoundError as e:
        st.error(str(e))
        return

    bp, cp = _edge_pairs(boot), _edge_pairs(cv)
    both, boot_only, cv_only = bp & cp, bp - cp, cp - bp

    c1, c2, c3 = st.columns(3)
    c1.metric("In both pipelines", len(both))
    c2.metric("Bootstrap only", len(boot_only))
    c3.metric("CV only", len(cv_only))

    total = len(bp | cp)
    if total:
        st.caption(f"Agreement: **{len(both) / total:.0%}** of all DAG edges appear in both.")

    t1, t2, t3 = st.tabs(["In both", "Bootstrap only", "CV only"])
    with t1:
        st.dataframe(_pretty(both), hide_index=True, use_container_width=True)
    with t2:
        st.dataframe(_pretty(boot_only), hide_index=True, use_container_width=True)
    with t3:
        st.dataframe(_pretty(cv_only), hide_index=True, use_container_width=True)
