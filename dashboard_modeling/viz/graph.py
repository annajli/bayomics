"""Pure transforms from edge/mapping DataFrames to graph structures.

filter_edges / root_aliases are Streamlit-free and fully unit-tested.
build_agraph returns streamlit-agraph Node/Edge objects for the Explorer page.
"""
import pandas as pd
from streamlit_agraph import Node, Edge

ROOT_COLOR = "#C44E52"
NODE_COLOR = "#4C72B0"
ROOT_SIZE = 26
NODE_SIZE = 14


def filter_edges(edges: pd.DataFrame, threshold: float, dag_only: bool) -> pd.DataFrame:
    """Return the subset of edges to render.

    dag_only=True  -> only edges in the final averaged DAG (ignores threshold).
    dag_only=False -> all edges with strength >= threshold.
    """
    if dag_only:
        return edges[edges["in_final_dag"]].copy()
    return edges[edges["strength"] >= threshold].copy()


def root_aliases(mapping: pd.DataFrame) -> set:
    """Aliases of the discrete root nodes (age_group, sex, cmv)."""
    return set(mapping.loc[mapping["type"] == "discrete_root", "alias"])


def build_agraph(edges: pd.DataFrame, roots: set):
    """Build (nodes, edges) lists for streamlit-agraph from a filtered edge frame."""
    node_ids = set(edges["from"]) | set(edges["to"])
    nodes = []
    for nid in sorted(node_ids):
        is_root = nid in roots
        nodes.append(Node(
            id=nid,
            label=nid,
            size=ROOT_SIZE if is_root else NODE_SIZE,
            color=ROOT_COLOR if is_root else NODE_COLOR,
        ))
    agedges = []
    for _, r in edges.iterrows():
        strength = float(r["strength"])
        agedges.append(Edge(
            source=str(r["from"]),
            target=str(r["to"]),
            width=1 + 4 * strength,  # thicker = more stable edge
            title=f"strength={strength:.3f}, direction={float(r['direction']):.2f}",
        ))
    return nodes, agedges
