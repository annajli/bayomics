"""Single source of truth mapping each network level to its canonical artifact files.

The R outputs carry version sprawl (L1a uses `_v3`, L2 uses `_v2`) and the two
pipelines use different naming (bn_learning keeps the version suffix; bn_cv does not).
Nothing else in the app hardcodes a path — every page reads from here.
"""
from pathlib import Path

# config/ -> dashboard_modeling/ -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]
BN = REPO_ROOT / "bn_learning" / "output"
CV = REPO_ROOT / "bn_cv" / "output"

_CSV = BN / "networks" / "csv_networks"
_MAP = BN / "mappings"
_PNG = BN / "networks" / "visualized_networks"
_CVCSV = CV / "networks" / "csv_networks"

LEVEL_ORDER = ["L1a", "L1b", "L2", "L3a", "L3b", "L4", "L5", "L6", "L_all"]

NETWORKS = {
    "L1a": {
        "label": "L1a — Clinical (full panel)",
        "n_nodes": 37,
        "edges": _CSV / "L1a_clinical_full_v3_literature_edges.csv",
        "mapping": _MAP / "L1a_clinical_full_v3_node_mapping.csv",
        "png": _PNG / "L1a_clinical_full_v3_avg_opt.png",
        "cv_edges": _CVCSV / "L1a_clinical_full_literature_edges_cv.csv",
        "notes": "",
    },
    "L1b": {
        "label": "L1b — Clinical (curated)",
        "n_nodes": 25,
        "edges": _CSV / "L1b_clinical_curated_literature_edges.csv",
        "mapping": _MAP / "L1b_clinical_curated_node_mapping.csv",
        "png": _PNG / "L1b_clinical_curated_avg_opt.png",
        "cv_edges": _CVCSV / "L1b_clinical_curated_literature_edges_cv.csv",
        "notes": "Cleanest network.",
    },
    "L2": {
        "label": "L2 — Olink proteomics",
        "n_nodes": 37,
        "edges": _CSV / "L2_olink_v2_literature_edges.csv",
        "mapping": _MAP / "L2_olink_v2_node_mapping.csv",
        "png": _PNG / "L2_olink_v2_avg_opt.png",
        "cv_edges": _CVCSV / "L2_olink_literature_edges_cv.csv",
        "notes": "",
    },
    "L3a": {
        "label": "L3a — Clinical (full) + Olink",
        "n_nodes": 71,
        "edges": _CSV / "L3a_clinical_full_olink_literature_edges.csv",
        "mapping": _MAP / "L3a_clinical_full_olink_node_mapping.csv",
        "png": _PNG / "L3a_clinical_full_olink_avg_opt.png",
        "cv_edges": _CVCSV / "L3a_clinical_full_olink_literature_edges_cv.csv",
        "notes": "Contains some derived-column noise (see spec §11).",
    },
    "L3b": {
        "label": "L3b — Clinical (curated) + Olink",
        "n_nodes": 59,
        "edges": _CSV / "L3b_clinical_curated_olink_literature_edges.csv",
        "mapping": _MAP / "L3b_clinical_curated_olink_node_mapping.csv",
        "png": _PNG / "L3b_clinical_curated_olink_avg_opt.png",
        "cv_edges": _CVCSV / "L3b_clinical_curated_olink_literature_edges_cv.csv",
        "notes": "Best combined network.",
    },
    "L4": {
        "label": "L4 — Whole-blood pathways",
        "n_nodes": 53,
        "edges": _CSV / "L4_wb_pathways_literature_edges.csv",
        "mapping": _MAP / "L4_wb_pathways_node_mapping.csv",
        "png": _PNG / "L4_wb_pathways_avg_opt.png",
        "cv_edges": _CVCSV / "L4_wb_pathways_literature_edges_cv.csv",
        "notes": "",
    },
    "L5": {
        "label": "L5 — Cell-frequency immunophenotype",
        "n_nodes": 91,
        "edges": _CSV / "L5_freq_immunophenotype_literature_edges.csv",
        "mapping": _MAP / "L5_freq_immunophenotype_node_mapping.csv",
        "png": _PNG / "L5_freq_immunophenotype_avg_opt.png",
        "cv_edges": _CVCSV / "L5_freq_immunophenotype_literature_edges_cv.csv",
        "notes": "",
    },
    "L6": {
        "label": "L6 — Pseudobulk signaling",
        "n_nodes": 147,
        "edges": _CSV / "L6_pb_signaling_literature_edges.csv",
        "mapping": _MAP / "L6_pb_signaling_node_mapping.csv",
        "png": _PNG / "L6_pb_signaling_avg_opt.png",
        "cv_edges": _CVCSV / "L6_pb_signaling_literature_edges_cv.csv",
        "notes": "Large — use DAG-only view or a high strength threshold.",
    },
    "L_all": {
        "label": "L_all — Kitchen sink (all modalities)",
        "n_nodes": 352,
        "edges": _CSV / "L_all_literature_edges.csv",
        "mapping": _MAP / "L_all_node_mapping.csv",
        "png": _PNG / "L_all_avg_opt.png",
        "cv_edges": _CVCSV / "L_all_literature_edges_cv.csv",
        "notes": "352 nodes / ~123k scored edge rows. Interactive view is DAG-only "
                 "by default; use the PNG toggle for the full picture.",
    },
}
