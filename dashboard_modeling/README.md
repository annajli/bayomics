# Bayomics — Modeling Results Dashboard

Standalone Streamlit app presenting the Bayesian-network modeling results
(network ladder L1a–L_all, root→modality findings, bootstrap-vs-CV comparison,
flu-response validation). Display-only: reads exported CSV/PNG artifacts from
`../bn_learning/output/` and `../bn_cv/output/`. No R at runtime.

Separate from the EDA dashboard in `../dashboard/`.

## Run locally (from repo root)

```bash
pip install -r dashboard_modeling/requirements.txt
streamlit run dashboard_modeling/app.py
```

## Tests (from repo root)

```bash
pip install -r dashboard_modeling/requirements-dev.txt
PYTHONPATH=dashboard_modeling pytest dashboard_modeling/tests -v
```

## Deploy (Streamlit Community Cloud)

- Main file path: `dashboard_modeling/app.py`
- Requirements: `dashboard_modeling/requirements.txt`
- All data is in-repo (~50 MB); no external store needed.
