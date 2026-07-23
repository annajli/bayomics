import streamlit as st
from streamlit_agraph import agraph, Config

from config.network_registry import NETWORKS, LEVEL_ORDER
from data.loaders import load_edges, load_mapping
from viz.graph import filter_edges, root_aliases, build_agraph
from viz.charts import strength_histogram


def render():
    st.title("🔗 Network Explorer")

    level = st.selectbox(
        "Network level", LEVEL_ORDER,
        format_func=lambda lvl: NETWORKS[lvl]["label"],
    )
    meta = NETWORKS[level]
    if meta["notes"]:
        st.caption(f"ℹ️ {meta['notes']}")

    edges = load_edges(level)
    mapping = load_mapping(level)
    roots = root_aliases(mapping)

    # --- controls -----------------------------------------------------------
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        dag_only = st.toggle("Final DAG only", value=True,
                             help="Show just the averaged-network edges (recommended).")
    with c2:
        threshold = st.slider("Min edge strength", 0.0, 1.0, 0.85, 0.01,
                              disabled=dag_only,
                              help="Only used when 'Final DAG only' is off.")
    with c3:
        show_png = st.toggle("Publication view (PNG)", value=False,
                             help="Swap the live graph for the pre-rendered Rgraphviz image.")

    if show_png:
        if meta["png"].exists():
            st.image(str(meta["png"]), use_container_width=True,
                     caption=f"{meta['label']} — averaged DAG (avg_opt)")
        else:
            st.error(f"PNG not found: {meta['png']}")
        return

    sub = filter_edges(edges, threshold=threshold, dag_only=dag_only)
    st.caption(f"Showing **{len(sub)}** edges across "
               f"**{len(set(sub['from']) | set(sub['to']))}** nodes "
               f"(roots highlighted in red).")

    # --- interactive graph --------------------------------------------------
    nodes, agedges = build_agraph(sub, roots)
    config = Config(width=900, height=600, directed=True,
                    physics=True, hierarchical=False,
                    nodeHighlightBehavior=True, highlightColor="#F2C14E",
                    collapsible=False)
    agraph(nodes=nodes, edges=agedges, config=config)

    # --- edge table + histogram --------------------------------------------
    left, right = st.columns([3, 2])
    with left:
        st.markdown("##### Edge table")
        st.dataframe(
            sub[["from", "to", "strength", "direction", "in_final_dag"]]
            .sort_values("strength", ascending=False),
            hide_index=True, use_container_width=True, height=360,
        )
    with right:
        st.markdown("##### Strength distribution (all scored edges)")
        st.plotly_chart(strength_histogram(edges), use_container_width=True)
