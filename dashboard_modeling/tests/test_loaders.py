from data.loaders import load_edges, load_mapping, load_validation


def test_load_edges_has_expected_columns():
    df = load_edges("L2")
    for col in ("from", "to", "strength", "direction", "in_final_dag"):
        assert col in df.columns
    assert df["strength"].between(0, 1).all()
    assert df["in_final_dag"].dtype == bool


def test_load_edges_cv_variant():
    df = load_edges("L2", pipeline="cv")
    assert "strength" in df.columns
    assert len(df) > 0


def test_load_mapping_flags_roots():
    m = load_mapping("L2")
    roots = m.loc[m["type"] == "discrete_root", "alias"].tolist()
    assert set(roots) == {"age_group", "sex", "cmv"}


def test_load_validation_shape():
    v = load_validation()
    assert set(v["target"]) == {"seroconvert_any", "igg_log2fc_mean"}
    assert v["agrees"].dtype == bool
