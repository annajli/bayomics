"""
SoundLife Capstone — EDA Dashboard
==================================
A Streamlit dashboard that brings together every piece of the exploratory
data analysis for the SoundLife healthy-cohort capstone:

  • Differential Expression (DESeq2)  -> fully interactive (raw CSVs loaded live)
  • Clinical Labs & Metadata          -> rendered notebook (narrative + plots)
  • Flu Serology & HAI                -> rendered notebook (narrative + plots)
  • Pseudobulk QC report              -> embedded HTML report
  • Plasma & Whole Blood              -> embedded HTML report

Run with:
    streamlit run app.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import plotly.express as px
import plotly.graph_objects as go

# --------------------------------------------------------------------------- #
# Paths & configuration
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
EDA_DIR = HERE.parent / "eda"                       # the original EDA folder
DESEQ_DIR = (
    HERE.parent.parent
    / "UVA_Data_Capstone" / "Datasets" / "Differential_Expression"
)

PADJ_SIG = 0.05      # significance threshold on adjusted p-value
LFC_THR = 1.0        # |log2FC| threshold for a "strong / biologically meaningful" hit

# Each contrast: (filename, scRNA column, background level, foreground level, color)
COMPARISONS = {
    "Age Group": dict(
        file="sound-life_age-group_deseq2-results_2025-02-07.csv",
        col="cohort.cohortGuid", bg="BR1 Younger Adult", fg="BR2 Older Adult",
        color="#4C72B0"),
    "Biological Sex": dict(
        file="sound-life_biological-sex_deseq2-results_2025-02-07.csv",
        col="subject.biologicalSex", bg="Female", fg="Male",
        color="#DD8452"),
    "CMV Status": dict(
        file="sound-life_cmv-status_deseq2-results_2025-02-07.csv",
        col="subject.cmv", bg="Negative", fg="Positive",
        color="#55A868"),
    "Flu Vaccine": dict(
        file="sound-life_flu-vaccine_deseq2-results_2025-02-07.csv",
        col="sample.visitName", bg="Flu Year 1 Day 0", fg="Flu Year 1 Day 7",
        color="#C44E52"),
}

st.set_page_config(
    page_title="SoundLife Capstone — EDA",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --------------------------------------------------------------------------- #
# Data loading (cached)
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False)
def load_comparison(name: str) -> pd.DataFrame:
    """Load and annotate one DESeq2 results file."""
    meta = COMPARISONS[name]
    df = pd.read_csv(DESEQ_DIR / meta["file"])
    for c in ("padj", "pvalue", "log2fc", "stat"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["sig"] = df["padj"] < PADJ_SIG
    df["sig_strong"] = df["sig"] & (df["log2fc"].abs() >= LFC_THR)
    df["direction"] = np.where(df["log2fc"] > 0, "Up in Foreground", "Up in Background")
    return df


@st.cache_data(show_spinner=False)
def load_all() -> dict:
    return {name: load_comparison(name) for name in COMPARISONS}


@st.cache_data(show_spinner=False)
def overview_table(_data: dict) -> pd.DataFrame:
    rows = []
    for name, df in _data.items():
        m = COMPARISONS[name]
        n = len(df)
        rows.append(dict(
            Comparison=name,
            Column=m["col"],
            Background=m["bg"],
            Foreground=m["fg"],
            Tests=n,
            CellTypes=df["AIFI_L3"].nunique(),
            Genes=df["gene"].nunique(),
            Significant=int(df["sig"].sum()),
            **{"% Sig": round(100 * df["sig"].sum() / n, 1)},
            Strong=int(df["sig_strong"].sum()),
            Up=int((df["sig"] & (df["log2fc"] > 0)).sum()),
            Down=int((df["sig"] & (df["log2fc"] < 0)).sum()),
        ))
    return pd.DataFrame(rows)


def deseq_available() -> bool:
    return DESEQ_DIR.exists() and any(DESEQ_DIR.glob("*.csv"))


def embed_html(path: Path, height: int = 900):
    """Embed a self-contained HTML report inside a scrollable iframe."""
    if not path.exists():
        st.error(f"Report not found: `{path}`")
        return
    html = path.read_text(encoding="utf-8", errors="ignore")
    components.html(html, height=height, scrolling=True)


# --------------------------------------------------------------------------- #
# Page: Home / Overview
# --------------------------------------------------------------------------- #
def page_home():
    st.title("🧬 SoundLife Capstone — Exploratory Data Analysis")
    st.markdown(
        "Interactive companion to the SoundLife healthy-cohort capstone EDA. "
        "The dataset profiles a healthy adult cohort across **single-cell RNA-seq**, "
        "**clinical labs**, and **flu-vaccine serology**, with the goal of "
        "characterizing the immune landscape of aging."
    )

    if deseq_available():
        data = load_all()
        ov = overview_table(data)
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Contrasts", len(COMPARISONS))
        c2.metric("Total DEG tests", f"{int(ov['Tests'].sum()):,}")
        c3.metric("Significant (padj<0.05)", f"{int(ov['Significant'].sum()):,}")
        c4.metric("Strong (& |log2FC|≥1)", f"{int(ov['Strong'].sum()):,}")
    else:
        st.warning(
            "DESeq2 CSVs were not found on this machine, so the Differential "
            "Expression page will fall back to the static HTML report. "
            f"Expected at: `{DESEQ_DIR}`"
        )

    st.subheader("What's in this dashboard")
    st.markdown(
        """
| Section | Source | Type |
|---|---|---|
| **Differential Expression** | 4 DESeq2 result tables | 🔵 Interactive |
| **Clinical Labs & Metadata** | `clinical_labs_metadata.ipynb` | 📄 Rendered notebook |
| **Flu Serology & HAI** | `flu_serology_hai.ipynb` | 📄 Rendered notebook |
| **Pseudobulk QC Report** | `EDA_pseudobulk_report.html` | 📄 Embedded report |
| **Plasma & Whole Blood** | `plasma_wholeblood.html` | 📄 Embedded report |

Use the sidebar to navigate. The Differential Expression page is fully
interactive — filter by contrast, cell type, and gene, and explore volcano
plots, top hits, and the cell-type signal landscape.
        """
    )

    with st.expander("About the contrasts (DESeq2)"):
        st.markdown(
            "Each test = one **gene** in one **cell type** (`AIFI_L3`) for one "
            "contrast. A positive `log2FC` means higher expression in the "
            "**foreground** group."
        )
        st.dataframe(
            pd.DataFrame([
                dict(Comparison=k, Column=v["col"],
                     Background=v["bg"], Foreground=v["fg"])
                for k, v in COMPARISONS.items()
            ]),
            hide_index=True, use_container_width=True,
        )


# --------------------------------------------------------------------------- #
# Page: Differential Expression (interactive)
# --------------------------------------------------------------------------- #
def page_diffexp():
    st.title("🔬 Differential Expression (DESeq2)")

    if not deseq_available():
        st.info("Raw DESeq2 CSVs not found — showing the pre-rendered report instead.")
        embed_html(EDA_DIR / "differential_expression_EDA.html", height=1000)
        return

    data = load_all()
    ov = overview_table(data)

    # ---- overview across all contrasts -----------------------------------
    st.subheader("Significant results by contrast")
    long = ov.melt(id_vars="Comparison", value_vars=["Significant", "Strong"],
                   var_name="Type", value_name="Count")
    fig = px.bar(long, x="Comparison", y="Count", color="Type", barmode="group",
                 color_discrete_map={"Significant": "#9aa7b3", "Strong": "#1f3a5f"},
                 text="Count")
    fig.update_traces(textposition="outside")
    fig.update_layout(height=380, margin=dict(t=30, b=10))
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("Dataset overview table", expanded=False):
        st.dataframe(ov, hide_index=True, use_container_width=True)

    st.divider()

    # ---- single-contrast explorer ----------------------------------------
    st.subheader("Explore a contrast")
    col_sel, col_ct = st.columns([1, 2])
    with col_sel:
        name = st.selectbox("Contrast", list(COMPARISONS.keys()))
    df = data[name]
    meta = COMPARISONS[name]

    with col_ct:
        cell_types = sorted(df["AIFI_L3"].dropna().unique())
        chosen_ct = st.multiselect(
            "Cell types (AIFI_L3) — empty = all",
            cell_types, default=[],
            help="Filter the volcano plot and tables to specific cell types.")

    sub = df if not chosen_ct else df[df["AIFI_L3"].isin(chosen_ct)]

    st.caption(
        f"**{meta['fg']}** (foreground) vs **{meta['bg']}** (background) "
        f"· column `{meta['col']}` · "
        f"{int(sub['sig'].sum()):,} significant of {len(sub):,} tests "
        f"({int(sub['sig_strong'].sum()):,} strong)."
    )

    # ---- volcano ----------------------------------------------------------
    plot_df = sub.dropna(subset=["padj", "log2fc"]).copy()
    plot_df = plot_df[plot_df["padj"] > 0]
    plot_df["neglog10padj"] = -np.log10(plot_df["padj"])

    def category(r):
        if not r["sig"]:
            return "Not significant"
        return "Up in Foreground" if r["log2fc"] > 0 else "Up in Background"
    plot_df["Category"] = plot_df.apply(category, axis=1)

    # subsample the non-significant cloud so the plot renders quickly
    ns = plot_df[plot_df["Category"] == "Not significant"]
    if len(ns) > 8000:
        ns = ns.sample(8000, random_state=0)
    sig_pts = plot_df[plot_df["Category"] != "Not significant"]
    vdf = pd.concat([ns, sig_pts], ignore_index=True)

    fig = px.scatter(
        vdf, x="log2fc", y="neglog10padj", color="Category",
        hover_data={"gene": True, "AIFI_L3": True, "log2fc": ":.2f",
                    "padj": ":.2e", "neglog10padj": False},
        color_discrete_map={
            "Not significant": "#d3d8de",
            "Up in Foreground": "#C44E52",
            "Up in Background": "#4C72B0"},
        opacity=0.6,
    )
    fig.add_hline(y=-np.log10(PADJ_SIG), line_dash="dash", line_color="#555")
    fig.add_vline(x=LFC_THR, line_dash="dot", line_color="#999")
    fig.add_vline(x=-LFC_THR, line_dash="dot", line_color="#999")
    fig.update_traces(marker=dict(size=5))
    fig.update_layout(height=520, xaxis_title="log2(Fold Change)",
                      yaxis_title="-log10(padj)", legend_title="",
                      margin=dict(t=20))
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Dashed line = padj 0.05 · dotted lines = ±1 log2FC. "
               "Hover a point for gene / cell type / stats.")

    st.divider()

    # ---- top hits + cell-type ranking ------------------------------------
    left, right = st.columns([3, 2])
    with left:
        st.markdown("##### Top significant hits")
        gene_q = st.text_input("Search gene (optional)", "").strip().upper()
        top = sub[sub["sig"]].dropna(subset=["padj"]).copy()
        if gene_q:
            top = top[top["gene"].str.upper().str.contains(gene_q, na=False)]
        top = top.sort_values(["padj", "pvalue"]).head(200)
        show = top[["AIFI_L3", "gene", "log2fc", "padj", "direction"]].rename(
            columns={"AIFI_L3": "Cell type", "gene": "Gene",
                     "log2fc": "log2FC", "padj": "padj",
                     "direction": "Direction"})
        st.dataframe(
            show.style.format({"log2FC": "{:.2f}", "padj": "{:.2e}"}),
            hide_index=True, use_container_width=True, height=420)

    with right:
        st.markdown("##### Most differentially active cell types")
        ct_counts = (sub[sub["sig_strong"]].groupby("AIFI_L3").size()
                     .sort_values(ascending=False).head(15)
                     .rename("Strong DEGs").reset_index()
                     .rename(columns={"AIFI_L3": "Cell type"}))
        if len(ct_counts):
            figc = px.bar(ct_counts.iloc[::-1], x="Strong DEGs", y="Cell type",
                          orientation="h", color_discrete_sequence=[meta["color"]])
            figc.update_layout(height=420, margin=dict(t=10, l=10))
            st.plotly_chart(figc, use_container_width=True)
        else:
            st.info("No strong DEGs for this selection.")


# --------------------------------------------------------------------------- #
# Pages: rendered notebooks & embedded reports
# --------------------------------------------------------------------------- #
def page_clinical():
    st.title("🩸 Clinical Labs & Metadata")
    st.caption("Rendered from `clinical_labs_metadata.ipynb` "
               "(narrative + figures; code hidden).")
    embed_html(ASSETS / "clinical_labs_metadata.html", height=1100)


def page_serology():
    st.title("💉 Flu Serology & HAI")
    st.caption("Rendered from `flu_serology_hai.ipynb` "
               "(narrative + figures; code hidden).")
    embed_html(ASSETS / "flu_serology_hai.html", height=1100)


def page_pseudobulk():
    st.title("🧫 Pseudobulk QC Report")
    st.caption("Embedded report: `EDA_pseudobulk_report.html`")
    embed_html(EDA_DIR / "EDA_pseudobulk_report.html", height=1100)


def page_plasma():
    st.title("🧪 Plasma & Whole Blood")
    st.caption("Embedded report: `plasma_wholeblood.html`")
    embed_html(EDA_DIR / "plasma_wholeblood.html", height=1100)


# --------------------------------------------------------------------------- #
# Navigation
# --------------------------------------------------------------------------- #
pages = [
    st.Page(page_home, title="Overview", icon="🏠", default=True),
    st.Page(page_diffexp, title="Differential Expression", icon="🔬"),
    st.Page(page_clinical, title="Clinical Labs & Metadata", icon="🩸"),
    st.Page(page_serology, title="Flu Serology & HAI", icon="💉"),
    st.Page(page_pseudobulk, title="Pseudobulk Report", icon="🧫"),
    st.Page(page_plasma, title="Plasma & Whole Blood", icon="🧪"),
]
st.navigation(pages).run()
