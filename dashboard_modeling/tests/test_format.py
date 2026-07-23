import pandas as pd

from viz.format import agreement_rate, scorecard_display

SCORE = pd.DataFrame({
    "target": ["seroconvert_any", "seroconvert_any", "igg_log2fc_mean", "igg_log2fc_mean"],
    "term": ["age_groupOlder Adult", "sexMale", "age_groupOlder Adult", "sexMale"],
    "estimate": [0.415, -0.269, -0.060, -0.102],
    "direction": ["higher", "lower", "lower", "lower"],
    "expected_dir": ["lower", "lower", "lower", "lower"],
    "agrees": [False, True, True, True],
    "p_value": [0.496, 0.665, 0.636, 0.412],
})


def test_agreement_rate():
    assert agreement_rate(SCORE) == 0.75  # 3 of 4 agree


def test_scorecard_display_adds_symbol_column():
    disp = scorecard_display(SCORE)
    assert "Agreement" in disp.columns
    symbols = disp["Agreement"].tolist()
    assert symbols[0] == "✗" and symbols[1] == "✓"
