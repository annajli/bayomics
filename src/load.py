import pandas as pd


def load_csv(path) -> pd.DataFrame:
    """Load a CSV with GUID columns forced to str to prevent numeric coercion."""
    return pd.read_csv(
        path,
        dtype={"subject.subjectGuid": str, "sample.sampleKitGuid": str},
    )
