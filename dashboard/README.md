# Sound Life Cohort — BN Foundation Model Dashboard

Interactive Streamlit dashboard for exploring Bayesian-network
structure-learning results for the Sound Life Cohort project (layers
L1a, L1b, L2, L3a, L3b, L4, L5, L6, and the combined L_all network).

Every layer was learned under **two independent validation strategies**
from the same candidate node/arc set:

- **Cross-validation (CV)** — k-fold cross-validated structure learning,
  averaged network at an algorithm-chosen optimal arc-strength threshold.
- **Bootstrap** — bootstrap-resampled structure learning, averaged
  network at a conservative arc-strength threshold (with the shared 0.50
  "optimal-fit" line shown for reference).

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).

## Pages

- **Overview** — project summary, success criteria, layer-level stats.
  A sidebar toggle switches the whole page between the CV and Bootstrap
  model.
- **Network Explorer** — the interactive part. Pick a layer and a model
  (CV or Bootstrap), then either view the official averaged final DAG,
  or uncheck that box and sweep an arc-strength threshold yourself. The
  graph is a real `streamlit-agraph` widget: drag nodes, zoom/pan, click
  a node (or use the search box) to see its parents/children and
  download the filtered edge list.
- **Diagnostics** — interactive, rebuilt-from-data arc-strength
  histogram per layer for the selected model, plus a static reference
  image (the original held-out CV loss histogram for CV, or the
  original bootstrap edge-strength histogram with threshold lines for
  Bootstrap).
- **CV vs Bootstrap** — a dedicated head-to-head comparison page:
  cross-layer agreement summary (shared edges, direction-agreement,
  CV-only / Bootstrap-only edges, Jaccard similarity), summary charts,
  and a per-layer deep-dive with side-by-side reference networks,
  side-by-side histograms, and tables of exactly which arcs each model
  agrees, disagrees on direction, or finds uniquely.
- **Model Comparison** — the original CV-pipeline algorithm comparison
  (HC / tabu / MM-HC / bootstrap-averaged thresholds) and CV loss / BIC
  across layers, from the original model-comparison report. This
  benchmarks algorithms *within* the CV run — for CV-vs-Bootstrap, see
  the dedicated page above.
- **About / Methods** — project methodology, assumptions, and an
  explicit note on what the current edge tables do and don't contain
  yet (see below).

## Data

- `data/edges/cv/*.csv` and `data/edges/bootstrap/*.csv` — per-layer
  arc tables (`from, to, strength, direction, in_final_dag`) for the CV
  and Bootstrap runs respectively. Both share the identical candidate
  node/arc set per layer, so they are directly comparable row-for-row.
- `data/images/network/cv/*.png` — original static CV network
  renderings (averaged + tabu).
- `data/images/network/bootstrap/*.png` — original static Bootstrap
  network renderings (averaged at the optimal threshold, averaged at
  the conservative threshold, and tabu).
- `data/images/hist/cv/*.png` — original CV edge-strength / held-out
  CV-loss histograms.
- `data/images/hist/bootstrap/*.png` — original Bootstrap edge-strength
  histograms, with the optimal-fit (0.50) and conservative threshold
  lines marked.
- `data/reports/model_comparison.csv` — CV-pipeline algorithm-comparison
  table, extracted from the original `03_model_comparison_cv.html`
  report.
- `data/reports/bootstrap_thresholds.csv` — per-layer optimal-fit
  (0.50) and conservative arc-strength thresholds used for the
  Bootstrap final DAGs, read from the original bootstrap histograms.

### Important note on "literature" edges

Despite some source filenames referencing "literature edges," these
files currently contain **structure-learning output** (arc strength,
direction probability, final-DAG membership) for CV and Bootstrap
respectively — not the per-edge literature support/novelty
classification (supported / plausible / novel / conflicting) described
in the project write-up. That classification is separate, still-pending
work. Once it exists, join it onto these tables on `(from, to)` and the
dashboard can add a literature-support filter/legend to Network
Explorer without any other structural changes.

### Layer name mapping (Bootstrap uploads → dashboard layer keys)

A few layers had multiple candidate uploads (versioned re-runs). The
version used in this dashboard was the one whose node set exactly
matched the corresponding CV layer's node set:

| Dashboard layer | Bootstrap source file prefix |
|---|---|
| L1a | `L1a_clinical_full_v3` |
| L1b | `L1b_clinical_curated` |
| L2  | `L2_olink_v2` |
| L3a | `L3a_clinical_full_olink` |
| L3b | `L3b_clinical_curated_olink` |
| L4  | `L4_wb_pathways` |
| L5  | `L5_freq_immunophenotype` |
| L6  | `L6_pb_signaling` |
| L_all | `L_all` |

## Deploying

This is a standard Streamlit app — deploy as-is to Streamlit Community
Cloud, or any host that can run `streamlit run app.py` with the packages
in `requirements.txt`. No secrets or external API calls are required; all
data ships in `data/`.
