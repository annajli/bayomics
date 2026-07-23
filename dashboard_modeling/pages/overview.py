import pandas as pd
import streamlit as st

from config.network_registry import NETWORKS, LEVEL_ORDER
from data.loaders import load_validation
from viz.format import agreement_rate


def render():
    st.title("🧬 Bayomics — Bayesian Network Modeling Results")
    st.markdown(
        "Interpretable joint model of a healthy adult immune system, learned from "
        "**multi-omics baseline data** (Sound Life cohort, n≈92). This dashboard presents "
        "the modeling phase: a ladder of Bayesian networks from single modalities up to a "
        "352-node kitchen-sink model, the robust **root→modality** biology, a "
        "bootstrap-vs-cross-validation robustness check, and hold-out **flu-response "
        "validation**."
    )

    ladder = pd.DataFrame([
        {"Level": lvl, "Network": NETWORKS[lvl]["label"], "Nodes": NETWORKS[lvl]["n_nodes"]}
        for lvl in LEVEL_ORDER
    ])
    try:
        val = load_validation()
        rate = agreement_rate(val)
    except FileNotFoundError:
        rate = None

    c1, c2, c3 = st.columns(3)
    c1.metric("Network levels", len(LEVEL_ORDER))
    c2.metric("Largest network", f"{NETWORKS['L_all']['n_nodes']} nodes")
    c3.metric("Validation agreement",
              "—" if rate is None else f"{rate:.0%}")

    st.subheader("The network ladder")
    st.caption("Each level is a Conditional Linear Gaussian BN over a modality slice; "
               "L_all combines everything. Explore any of them under **Network Explorer**.")
    st.dataframe(ladder, hide_index=True, use_container_width=True)

    st.subheader("How to read this dashboard")
    st.markdown(
        "- **Network Explorer** — interactive graphs from the bootstrap edge lists.\n"
        "- **Root → Modality Findings** — the headline biology (CMV, Age, Sex effects).\n"
        "- **Bootstrap vs CV** — edge stability vs predictive generalization.\n"
        "- **Validation** — does baseline immune architecture predict flu response?"
    )
