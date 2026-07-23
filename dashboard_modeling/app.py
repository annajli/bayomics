"""Bayomics modeling-results dashboard entry point.

Run from repo root:  streamlit run dashboard_modeling/app.py
On Streamlit Cloud, set the main file to dashboard_modeling/app.py.
"""
import sys
from pathlib import Path

# Make `config`, `data`, `viz`, `pages` importable when Streamlit runs this file
# directly (its own directory is not automatically on sys.path).
sys.path.insert(0, str(Path(__file__).resolve().parent))

import streamlit as st  # noqa: E402

from pages import (  # noqa: E402
    overview, network_explorer, root_findings, bootstrap_vs_cv, validation,
)

st.set_page_config(
    page_title="Bayomics — Modeling Results",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

nav = st.navigation([
    st.Page(overview.render, title="Overview", icon="🏠", default=True),
    st.Page(network_explorer.render, title="Network Explorer", icon="🔗"),
    st.Page(root_findings.render, title="Root → Modality Findings", icon="🧭"),
    st.Page(bootstrap_vs_cv.render, title="Bootstrap vs CV", icon="⚖️"),
    st.Page(validation.render, title="Validation", icon="✅"),
])
nav.run()
