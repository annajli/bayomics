import pandas as pd

from viz.graph import filter_edges, root_aliases, build_agraph

EDGES = pd.DataFrame({
    "from": ["age_group", "CXCL9", "IL6", "TNF"],
    "to":   ["GDF15", "CXCL10", "CRP", "HAVCR2"],
    "strength": [0.99, 0.90, 0.40, 0.88],
    "direction": [1.0, 0.52, 0.5, 0.61],
    "in_final_dag": [True, True, False, True],
})

MAPPING = pd.DataFrame({
    "alias": ["age_group", "sex", "cmv", "GDF15", "CXCL9"],
    "csv_column": ["subject.ageGroup", "subject.biologicalSex",
                   "cmv.igg_serology_interpretation", "olink.GDF15", "olink.CXCL9"],
    "type": ["discrete_root", "discrete_root", "discrete_root", "continuous", "continuous"],
    "log_transformed": [False, False, False, False, False],
})


def test_filter_edges_dag_only():
    out = filter_edges(EDGES, threshold=0.0, dag_only=True)
    assert len(out) == 3  # the three in_final_dag==True rows
    assert out["in_final_dag"].all()


def test_filter_edges_threshold():
    out = filter_edges(EDGES, threshold=0.89, dag_only=False)
    assert set(zip(out["from"], out["to"])) == {("age_group", "GDF15"), ("CXCL9", "CXCL10")}


def test_root_aliases():
    assert root_aliases(MAPPING) == {"age_group", "sex", "cmv"}


def test_build_agraph_counts_and_root_flagging():
    sub = filter_edges(EDGES, threshold=0.0, dag_only=True)
    nodes, edges = build_agraph(sub, roots={"age_group", "sex", "cmv"})
    node_ids = {n.id for n in nodes}
    assert node_ids == {"age_group", "GDF15", "CXCL9", "CXCL10", "TNF", "HAVCR2"}
    assert len(edges) == 3
    root_node = next(n for n in nodes if n.id == "age_group")
    non_root = next(n for n in nodes if n.id == "GDF15")
    assert root_node.size > non_root.size  # roots rendered larger
