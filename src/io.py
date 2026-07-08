import os
import pandas as pd

REQUIRED_INDEX_COLS = ["subject.subjectGuid", "sample.sampleKitGuid"]
MANIFEST_COLS = ["feature_name", "modality", "source_column", "rationale"]


def enforce_export_contract(df: pd.DataFrame, prefix: str) -> None:
    """Assert the output DataFrame meets the shared export contract."""
    for col in REQUIRED_INDEX_COLS:
        if col not in df.columns:
            raise ValueError(f"Missing required index column: {col}")
    feature_cols = [c for c in df.columns if c not in REQUIRED_INDEX_COLS]
    bad = [c for c in feature_cols if not c.startswith(prefix)]
    if bad:
        raise ValueError(f"Feature columns missing prefix '{prefix}': {bad}")


def write_processed(df: pd.DataFrame, modality_name: str, output_dir: str = "data/processed") -> str:
    """Write the wide output CSV to data/processed/ and return the path."""
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{modality_name}_baseline_wide.csv")
    df.to_csv(path, index=False)
    return path


def append_manifest(entries: list, manifest_path: str) -> None:
    """Append feature manifest entries, creating the file if it doesn't exist."""
    new_rows = pd.DataFrame(entries, columns=MANIFEST_COLS)
    if os.path.exists(manifest_path):
        existing = pd.read_csv(manifest_path)
        combined = pd.concat([existing, new_rows], ignore_index=True)
    else:
        combined = new_rows
    combined.to_csv(manifest_path, index=False)
