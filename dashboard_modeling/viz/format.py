"""Formatting helpers for the flu-response validation scorecard."""
import pandas as pd


def agreement_rate(score: pd.DataFrame) -> float:
    """Fraction of rows where observed direction agrees with the expected direction."""
    if len(score) == 0:
        return 0.0
    return round(score["agrees"].mean(), 4)


def scorecard_display(score: pd.DataFrame) -> pd.DataFrame:
    """Human-readable version of the scorecard for st.dataframe."""
    disp = score.copy()
    disp["Agreement"] = disp["agrees"].map({True: "✓", False: "✗"})
    disp = disp.rename(columns={
        "target": "Target",
        "term": "Root term",
        "estimate": "Estimate",
        "direction": "Observed",
        "expected_dir": "Expected",
        "p_value": "p-value",
    })
    return disp[["Target", "Root term", "Estimate", "Observed", "Expected",
                 "Agreement", "p-value"]]
