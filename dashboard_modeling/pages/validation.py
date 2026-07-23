import pandas as pd
import streamlit as st

from config.network_registry import NETWORKS
from data.loaders import load_edges, load_validation
from viz.format import agreement_rate, scorecard_display


def _has_literature_columns() -> bool:
    """True once the lit-review module has written verdict columns onto the edges.

    Reads only the CSV header (nrows=0) so it never pulls the full ~123k-row
    L_all edge list just to inspect column names.
    """
    path = NETWORKS["L_all"]["edges"]
    if not path.exists():
        return False
    cols = pd.read_csv(path, nrows=0).columns
    return any(c.lower().startswith("lit_") for c in cols)


def render():
    st.title("✅ Flu-Response Validation")
    st.markdown(
        "Does **baseline immune architecture predict post-vaccination flu response?** "
        "The fitted model's implied root effects on flu-response targets are compared "
        "against the expected direction from the literature, on held-out serology/HAI data."
    )

    try:
        score = load_validation()
    except FileNotFoundError as e:
        st.error(str(e))
        return

    st.metric("Direction agreement", f"{agreement_rate(score):.0%}",
              help="Share of target×root tests whose observed direction matches expectation.")

    st.subheader("Empirical scorecard")
    st.caption("Targets: `seroconvert_any` (binary response) and `igg_log2fc_mean` "
               "(continuous fold-change). ✓ = observed direction matches the literature.")
    disp = scorecard_display(score)
    st.dataframe(
        disp.style.format({"Estimate": "{:.3f}", "p-value": "{:.3f}"}, na_rep="—"),
        hide_index=True, use_container_width=True,
    )
    st.caption("Note: this is a small, mostly ceiling-limited cohort (overall seroconversion "
               "≈8%), so p-values are large — the value is in *direction agreement*, not "
               "significance.")

    st.divider()
    st.subheader("Literature-review-module verdicts")
    if _has_literature_columns():
        lit = load_edges("L_all")
        lit_cols = ["from", "to", "strength"] + [c for c in lit.columns
                                                 if c.lower().startswith("lit_")]
        st.dataframe(lit[lit_cols], hide_index=True, use_container_width=True)
    else:
        st.info("⏳ **Pending.** The AI-augmented literature-review module has not yet written "
                "edge verdicts onto the network edge lists. This section will populate "
                "automatically once `lit_*` columns appear on the edge CSVs.")
