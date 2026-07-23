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

# Every page function is named `render`, so Streamlit would infer the same URL
# pathname for all of them and error. Give each an explicit, unique url_path.
nav = st.navigation([
    st.Page(overview.render, title="Overview", icon="🏠", url_path="overview", default=True),
    st.Page(network_explorer.render, title="Network Explorer", icon="🔗", url_path="network-explorer"),
    st.Page(root_findings.render, title="Root → Modality Findings", icon="🧭", url_path="root-findings"),
    st.Page(bootstrap_vs_cv.render, title="Bootstrap vs CV", icon="⚖️", url_path="bootstrap-vs-cv"),
    st.Page(validation.render, title="Validation", icon="✅", url_path="validation"),
])
nav.run()
