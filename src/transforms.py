import re


def to_snake_case(name: str) -> str:
    """Convert a cell type name (may contain spaces/mixed case) to snake_case."""
    return re.sub(r"[^a-zA-Z0-9]+", "_", name).strip("_").lower()


# Stubs for later modalities
# def clr(counts: pd.DataFrame) -> pd.DataFrame: ...
# def ssgsea(adata, gene_sets: dict) -> pd.DataFrame: ...
