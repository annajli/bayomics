import pandas as pd


def assert_no_duplicates(df: pd.DataFrame, key_cols: list) -> None:
    """Raise ValueError if any combination of key_cols appears more than once."""
    dupes = df.duplicated(subset=key_cols, keep=False)
    if dupes.any():
        raise ValueError(
            f"Duplicate rows found for keys {key_cols}:\n"
            f"{df[dupes][key_cols].drop_duplicates().to_string()}"
        )


def flag_missingness(df: pd.DataFrame, threshold: float = 0.2) -> pd.DataFrame:
    """Return per-feature missingness summary; features > threshold flagged as exclude_candidate."""
    frac = df.isnull().mean()
    return pd.DataFrame(
        {
            "feature": frac.index,
            "missing_frac": frac.values,
            "exclude_candidate": frac.values > threshold,
        }
    ).reset_index(drop=True)
