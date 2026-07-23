"""Cached, Streamlit-aware readers for the modeling artifacts.

Each reader returns a clean DataFrame with normalized dtypes. Paths come only from
the network registry. When Streamlit is not running (e.g. pytest) the cache
decorator degrades to a no-op passthrough.
"""
import pandas as pd

from config.network_registry import NETWORKS, BN

try:
    import streamlit as st
    cache = st.cache_data(show_spinner=False)
except Exception:  # pragma: no cover - only hit outside a Streamlit runtime
    def cache(func):
        return func


def _to_bool(series: pd.Series) -> pd.Series:
    return (
        series.astype(str).str.strip().str.upper().map({"TRUE": True, "FALSE": False})
    )


@cache
def load_edges(level: str, pipeline: str = "bootstrap") -> pd.DataFrame:
    """Load one network's edge list. pipeline in {"bootstrap", "cv"}."""
    key = "edges" if pipeline == "bootstrap" else "cv_edges"
    path = NETWORKS[level][key]
    if not path.exists():
        raise FileNotFoundError(f"Edge list not found for {level} ({pipeline}): {path}")
    df = pd.read_csv(path)
    df["strength"] = pd.to_numeric(df["strength"], errors="coerce")
    df["direction"] = pd.to_numeric(df["direction"], errors="coerce")
    df["in_final_dag"] = _to_bool(df["in_final_dag"]).fillna(False).astype(bool)
    return df


@cache
def load_mapping(level: str) -> pd.DataFrame:
    """Load one network's node mapping (alias, csv_column, type, log_transformed)."""
    path = NETWORKS[level]["mapping"]
    if not path.exists():
        raise FileNotFoundError(f"Node mapping not found for {level}: {path}")
    return pd.read_csv(path)


@cache
def load_validation() -> pd.DataFrame:
    """Load the flu-response empirical/literature validation scorecard."""
    path = BN / "validation" / "emp_literature_scorecard.csv"
    if not path.exists():
        raise FileNotFoundError(f"Validation scorecard not found: {path}")
    df = pd.read_csv(path)
    df["agrees"] = _to_bool(df["agrees"]).fillna(False).astype(bool)
    df["estimate"] = pd.to_numeric(df["estimate"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    return df
