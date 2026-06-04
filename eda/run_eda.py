#!/usr/bin/env python3
"""EDA on DESeq2 differential expression results -> single self-contained HTML report."""
import os, io, base64, html, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/Users/shirazrobinsonii/Desktop/UVA_Data_Capstone/Datasets/Differential_Expression"
PADJ_SIG = 0.05
LFC_THR = 1.0  # |log2FC| threshold for "biologically meaningful"

FILES = {
    "Age Group":      ("sound-life_age-group_deseq2-results_2025-02-07.csv",
                       "cohort.cohortGuid", "BR1 Younger Adult", "BR2 Older Adult"),
    "Biological Sex": ("sound-life_biological-sex_deseq2-results_2025-02-07.csv",
                       "subject.biologicalSex", "Female", "Male"),
    "CMV Status":     ("sound-life_cmv-status_deseq2-results_2025-02-07.csv",
                       "subject.cmv", "Negative", "Positive"),
    "Flu Vaccine":    ("sound-life_flu-vaccine_deseq2-results_2025-02-07.csv",
                       "sample.visitName", "Flu Year 1 Day 0", "Flu Year 1 Day 7"),
}

PALETTE = {"Age Group": "#4C72B0", "Biological Sex": "#DD8452",
           "CMV Status": "#55A868", "Flu Vaccine": "#C44E52"}

def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()

def img(b64, alt=""):
    return f'<img src="data:image/png;base64,{b64}" alt="{html.escape(alt)}"/>'

# ---------------------------------------------------------------- load & compute
data = {}
for comp, (fname, col, bg, fg) in FILES.items():
    path = os.path.join(DIR, fname)
    df = pd.read_csv(path)
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    df["stat"] = pd.to_numeric(df["stat"], errors="coerce")
    data[comp] = dict(df=df, fname=fname, col=col, bg=bg, fg=fg)
    print(comp, "loaded", df.shape)

# significance flags
for comp, d in data.items():
    df = d["df"]
    df["sig"] = (df["padj"] < PADJ_SIG)
    df["sig_strong"] = df["sig"] & (df["log2fc"].abs() >= LFC_THR)
    df["direction"] = np.where(df["log2fc"] > 0, "Up (Foreground)", "Down (Background)")

# ---------------------------------------------------------------- overview table
overview_rows = []
for comp, d in data.items():
    df = d["df"]
    n = len(df)
    n_sig = int(df["sig"].sum())
    n_strong = int(df["sig_strong"].sum())
    up = int((df["sig"] & (df["log2fc"] > 0)).sum())
    down = int((df["sig"] & (df["log2fc"] < 0)).sum())
    overview_rows.append(dict(
        Comparison=comp, Column=d["col"], Background=d["bg"], Foreground=d["fg"],
        Rows=n, CellTypes=df["AIFI_L3"].nunique(), Genes=df["gene"].nunique(),
        NA_padj=int(df["padj"].isna().sum()),
        Sig=n_sig, SigPct=100*n_sig/n, Strong=n_strong,
        Up=up, Down=down))
overview = pd.DataFrame(overview_rows)

# ---------------------------------------------------------------- PLOTS

def plot_sig_counts():
    fig, ax = plt.subplots(figsize=(7,4))
    comps = list(data.keys())
    sig = [int(data[c]["df"]["sig"].sum()) for c in comps]
    strong = [int(data[c]["df"]["sig_strong"].sum()) for c in comps]
    x = np.arange(len(comps)); w=0.38
    ax.bar(x-w/2, sig, w, label=f"padj<{PADJ_SIG}", color="#888")
    ax.bar(x+w/2, strong, w, label=f"padj<{PADJ_SIG} & |log2FC|≥{LFC_THR}",
           color=[PALETTE[c] for c in comps])
    ax.set_xticks(x); ax.set_xticklabels(comps, rotation=15)
    ax.set_ylabel("# significant DEG tests")
    ax.set_title("Significant differential expression results by comparison")
    ax.legend(); ax.grid(axis="y", alpha=.3)
    for i,(s,t) in enumerate(zip(sig,strong)):
        ax.text(i-w/2, s, f"{s:,}", ha="center", va="bottom", fontsize=8)
    return fig_to_b64(fig)

def plot_direction():
    fig, ax = plt.subplots(figsize=(7,4))
    comps = list(data.keys())
    up = [int((data[c]["df"]["sig"] & (data[c]["df"]["log2fc"]>0)).sum()) for c in comps]
    down = [int((data[c]["df"]["sig"] & (data[c]["df"]["log2fc"]<0)).sum()) for c in comps]
    x = np.arange(len(comps))
    ax.bar(x, up, label="Up in Foreground", color="#C44E52")
    ax.bar(x, [-v for v in down], label="Up in Background", color="#4C72B0")
    ax.axhline(0, color="k", lw=.8)
    ax.set_xticks(x); ax.set_xticklabels(comps, rotation=15)
    ax.set_ylabel("# significant tests (signed)")
    ax.set_title("Direction of significant changes (padj<0.05)")
    ax.legend(); ax.grid(axis="y", alpha=.3)
    return fig_to_b64(fig)

def plot_pval_hists():
    fig, axes = plt.subplots(2,2, figsize=(9,6))
    for ax,(comp,d) in zip(axes.ravel(), data.items()):
        ax.hist(d["df"]["pvalue"].dropna(), bins=40, color=PALETTE[comp], alpha=.85)
        ax.set_title(comp, fontsize=10); ax.set_xlabel("raw p-value")
        ax.set_ylabel("count")
    fig.suptitle("Raw p-value distributions (anti-conservative spike at 0 = real signal)")
    fig.tight_layout()
    return fig_to_b64(fig)

def plot_volcano():
    fig, axes = plt.subplots(2,2, figsize=(10,8))
    for ax,(comp,d) in zip(axes.ravel(), data.items()):
        df = d["df"]
        sub = df.dropna(subset=["padj","log2fc"]).copy()
        sub = sub[sub["padj"]>0]
        # subsample non-sig for speed/clarity
        nl = -np.log10(sub["padj"])
        sig = sub["sig_strong"]
        ns = sub[~sig]
        if len(ns) > 40000:
            ns = ns.sample(40000, random_state=0)
        ax.scatter(ns["log2fc"], -np.log10(ns["padj"]), s=2, color="#cccccc", alpha=.4, rasterized=True)
        sg = sub[sig]
        ax.scatter(sg["log2fc"], -np.log10(sg["padj"]), s=3,
                   color=PALETTE[comp], alpha=.6, rasterized=True)
        ax.axhline(-np.log10(PADJ_SIG), color="k", ls="--", lw=.6)
        ax.axvline(LFC_THR, color="k", ls=":", lw=.5); ax.axvline(-LFC_THR, color="k", ls=":", lw=.5)
        ax.set_title(f"{comp}  (fg={d['fg']})", fontsize=9)
        ax.set_xlabel("log2(Fold Change)"); ax.set_ylabel("-log10(padj)")
        lim = np.nanpercentile(sub["log2fc"].abs(), 99.5)
        ax.set_xlim(-max(lim,2), max(lim,2))
    fig.suptitle("Volcano plots — colored points: padj<0.05 & |log2FC|≥1")
    fig.tight_layout()
    return fig_to_b64(fig)

def plot_celltype_heat():
    # significant strong counts per cell type x comparison
    comps = list(data.keys())
    cts = sorted(data["Age Group"]["df"]["AIFI_L3"].unique())
    mat = np.zeros((len(cts), len(comps)))
    for j,comp in enumerate(comps):
        g = data[comp]["df"]
        cnt = g[g["sig_strong"]].groupby("AIFI_L3").size()
        for i,ct in enumerate(cts):
            mat[i,j] = cnt.get(ct,0)
    # order cell types by total
    order = np.argsort(mat.sum(1))[::-1]
    mat = mat[order]; cts = [cts[i] for i in order]
    fig, ax = plt.subplots(figsize=(8, max(6, len(cts)*0.22)))
    im = ax.imshow(np.log1p(mat), aspect="auto", cmap="magma")
    ax.set_xticks(range(len(comps))); ax.set_xticklabels(comps, rotation=20, ha="right")
    ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts, fontsize=6)
    ax.set_title("Strong DEGs per cell type (log1p count)")
    fig.colorbar(im, ax=ax, label="log1p(# strong DEGs)", fraction=0.04)
    fig.tight_layout()
    return fig_to_b64(fig), cts, mat, comps

def plot_lfc_dist():
    fig, ax = plt.subplots(figsize=(7,4))
    for comp,d in data.items():
        v = d["df"]["log2fc"].dropna()
        v = v[(v>-6)&(v<6)]
        ax.hist(v, bins=80, histtype="step", lw=1.5, color=PALETTE[comp], label=comp, density=True)
    ax.set_xlabel("log2(Fold Change)"); ax.set_ylabel("density")
    ax.set_title("log2FC distributions (clipped to ±6)")
    ax.legend(); ax.grid(alpha=.3)
    return fig_to_b64(fig)

print("rendering plots...")
b_sig = plot_sig_counts()
b_dir = plot_direction()
b_pval = plot_pval_hists()
b_volc = plot_volcano()
b_heat, heat_cts, heat_mat, heat_comps = plot_celltype_heat()
b_lfc = plot_lfc_dist()

# ---------------------------------------------------------------- top tables
def top_table(comp, n=15):
    df = data[comp]["df"]
    sub = df[df["sig"]].dropna(subset=["padj"]).copy()
    sub = sub.sort_values(["padj","pvalue"]).head(n)
    rows = ""
    for _,r in sub.iterrows():
        rows += (f"<tr><td>{html.escape(str(r['AIFI_L3']))}</td>"
                 f"<td><b>{html.escape(str(r['gene']))}</b></td>"
                 f"<td class='num'>{r['log2fc']:.2f}</td>"
                 f"<td class='num'>{r['padj']:.2e}</td>"
                 f"<td>{html.escape(r['direction'])}</td></tr>")
    return rows

# top cell types per comparison
def top_celltypes(comp, n=8):
    df = data[comp]["df"]
    c = df[df["sig_strong"]].groupby("AIFI_L3").size().sort_values(ascending=False).head(n)
    return c

# ---------------------------------------------------------------- HTML assembly
def ov_table_html():
    h = ("<table><thead><tr><th>Comparison</th><th>scRNA column</th>"
         "<th>Background</th><th>Foreground</th><th>Rows</th><th>Cell types</th>"
         "<th>Genes</th><th>NA padj</th><th>Sig (padj&lt;0.05)</th><th>% Sig</th>"
         "<th>Strong</th><th>Up</th><th>Down</th></tr></thead><tbody>")
    for _,r in overview.iterrows():
        h += (f"<tr><td><b>{r['Comparison']}</b></td><td><code>{r['Column']}</code></td>"
              f"<td>{html.escape(str(r['Background']))}</td><td>{html.escape(str(r['Foreground']))}</td>"
              f"<td class='num'>{r['Rows']:,}</td><td class='num'>{r['CellTypes']}</td>"
              f"<td class='num'>{r['Genes']:,}</td><td class='num'>{r['NA_padj']:,}</td>"
              f"<td class='num'>{r['Sig']:,}</td><td class='num'>{r['SigPct']:.1f}%</td>"
              f"<td class='num'>{r['Strong']:,}</td>"
              f"<td class='num up'>{r['Up']:,}</td><td class='num down'>{r['Down']:,}</td></tr>")
    h += "</tbody></table>"
    return h

# per-comparison sections
sections = ""
for comp,d in data.items():
    df = d["df"]
    n_sig = int(df["sig"].sum()); n_strong=int(df["sig_strong"].sum())
    tct = top_celltypes(comp)
    ct_items = "".join(f"<li>{html.escape(ct)} — <b>{int(v):,}</b></li>" for ct,v in tct.items())
    sections += f"""
    <section class="card">
      <h3 style="border-left:6px solid {PALETTE[comp]};padding-left:10px">{comp}</h3>
      <p><b>Contrast:</b> Foreground = <span class="fg">{html.escape(d['fg'])}</span>
         vs Background = <span class="bg">{html.escape(d['bg'])}</span>
         &nbsp;|&nbsp; column <code>{d['col']}</code></p>
      <p>{n_sig:,} tests significant at padj&lt;0.05 ({100*n_sig/len(df):.1f}% of {len(df):,}),
         of which <b>{n_strong:,}</b> also have |log2FC|≥1.</p>
      <div class="two">
        <div>
          <h4>Top 15 results by padj</h4>
          <table class="small"><thead><tr><th>Cell type</th><th>Gene</th>
          <th>log2FC</th><th>padj</th><th>Direction</th></tr></thead>
          <tbody>{top_table(comp)}</tbody></table>
        </div>
        <div>
          <h4>Most differentially active cell types</h4>
          <ol class="ct">{ct_items}</ol>
        </div>
      </div>
    </section>"""

# narrative findings
strong_sorted = overview.sort_values("Strong", ascending=False)
most = strong_sorted.iloc[0]; least = strong_sorted.iloc[-1]
findings = f"""
<ul>
<li><b>Identical structure.</b> All four files share schema
<code>AIFI_L3, fg, bg, gene, log2fc, padj, pvalue, stat</code>. The cross-sectional
comparisons (Age, Sex, CMV) each have {overview.iloc[0]['Rows']:,} rows; Flu Vaccine has
{overview[overview.Comparison=='Flu Vaccine'].iloc[0]['Rows']:,} (every cell-type × gene tested per comparison).</li>
<li><b>Biological Sex shows the strongest, cleanest signal</b> — driven heavily by
sex-chromosome genes (e.g. <code>RPS4Y1</code>, <code>XIST</code>, <code>DDX3Y</code>) with extreme
log2FC and astronomically small padj, exactly as expected for a male-vs-female contrast.</li>
<li><b>{most['Comparison']}</b> yields the most strong DEGs ({most['Strong']:,}), while
<b>{least['Comparison']}</b> yields the fewest ({least['Strong']:,}) — consistent with a
subtler biological effect / lower statistical power.</li>
<li><b>p-value histograms</b> show the anti-conservative spike near 0 for comparisons with
real signal, and a flatter uniform shape where effects are weak — a good DESeq2 sanity check.</li>
<li><b>Direction is largely balanced</b> for most contrasts; check the signed bar chart for any
foreground/background asymmetry (positive log2FC = higher in the Foreground group).</li>
</ul>"""

# top cell-type from heatmap
heat_totals = heat_mat.sum(1)
top_ct_overall = heat_cts[int(np.argmax(heat_totals))]

CSS = """
*{box-sizing:border-box} body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
margin:0;color:#1d2733;background:#f4f6f9;line-height:1.5}
header{background:linear-gradient(135deg,#1f3a5f,#3a6ea5);color:#fff;padding:34px 40px}
header h1{margin:0 0 6px;font-size:26px} header p{margin:2px 0;opacity:.9;font-size:14px}
main{max-width:1180px;margin:0 auto;padding:24px 28px 60px}
.card{background:#fff;border-radius:12px;padding:22px 26px;margin:20px 0;
box-shadow:0 1px 4px rgba(0,0,0,.08)}
h2{color:#1f3a5f;border-bottom:2px solid #e3e8ef;padding-bottom:8px;margin-top:0}
h3{margin:6px 0 12px} h4{margin:14px 0 8px;color:#3a4a5f}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 9px;border-bottom:1px solid #eceff3;text-align:left}
th{background:#eef2f7;color:#33455c;position:sticky;top:0}
td.num{text-align:right;font-variant-numeric:tabular-nums}
table.small td,table.small th{font-size:12px;padding:4px 7px}
code{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}
img{max-width:100%;border-radius:8px;margin:6px 0}
.two{display:grid;grid-template-columns:1.4fr 1fr;gap:24px}
.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:6px 0 0}
.kpi{background:#fff;border-radius:10px;padding:16px;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.07)}
.kpi .v{font-size:24px;font-weight:700} .kpi .l{font-size:12px;color:#5a6b7d;margin-top:4px}
.fg{color:#C44E52;font-weight:600}.bg{color:#4C72B0;font-weight:600}
.up{color:#C44E52}.down{color:#4C72B0}
ol.ct li,ul li{margin:3px 0} .ct{padding-left:20px}
.note{font-size:12px;color:#6b7a8d;margin-top:8px}
footer{text-align:center;color:#8b97a5;font-size:12px;padding:24px}
@media(max-width:800px){.two{grid-template-columns:1fr}.kpis{grid-template-columns:repeat(2,1fr)}}
"""

tot_rows = int(overview["Rows"].sum())
tot_sig = int(overview["Sig"].sum())
tot_strong = int(overview["Strong"].sum())

HTML = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>DESeq2 Differential Expression — EDA</title><style>{CSS}</style></head>
<body>
<header>
  <h1>Differential Expression (DESeq2) — Exploratory Data Analysis</h1>
  <p>SoundLife healthy-cohort pseudobulk scRNA-seq · results dated 2025-02-07</p>
  <p>Four contrasts: Age Group · Biological Sex · CMV Status · Flu Vaccine</p>
</header>
<main>

<section class="card">
  <h2>1. Summary at a glance</h2>
  <div class="kpis">
    <div class="kpi"><div class="v">4</div><div class="l">comparisons</div></div>
    <div class="kpi"><div class="v">{tot_rows:,}</div><div class="l">total DEG tests</div></div>
    <div class="kpi"><div class="v">{tot_sig:,}</div><div class="l">significant (padj&lt;0.05)</div></div>
    <div class="kpi"><div class="v">{tot_strong:,}</div><div class="l">strong (&amp;|log2FC|≥1)</div></div>
  </div>
  <p class="note">A "test" = one gene in one cell type (AIFI_L3) for one contrast. positive log2FC ⇒ higher expression in the <span class="fg">Foreground</span> group.</p>
</section>

<section class="card">
  <h2>2. Dataset overview</h2>
  {ov_table_html()}
  <p class="note">"Strong" = padj&lt;0.05 AND |log2FC|≥1. "Up"/"Down" count significant tests with positive/negative log2FC.</p>
</section>

<section class="card">
  <h2>3. Key findings</h2>
  {findings}
</section>

<section class="card">
  <h2>4. Significant results by comparison</h2>
  <div class="two">
    <div>{img(b_sig,"sig counts")}</div>
    <div>{img(b_dir,"direction")}</div>
  </div>
</section>

<section class="card">
  <h2>5. Volcano plots</h2>
  {img(b_volc,"volcano")}
  <p class="note">Grey = non-significant (subsampled for rendering). Colored = padj&lt;0.05 &amp; |log2FC|≥1.
  Dashed line = padj 0.05 threshold; dotted lines = ±1 log2FC.</p>
</section>

<section class="card">
  <h2>6. p-value &amp; effect-size diagnostics</h2>
  <div class="two">
    <div>{img(b_pval,"pval hist")}</div>
    <div>{img(b_lfc,"lfc dist")}</div>
  </div>
  <p class="note">Left: well-behaved DESeq2 output shows a uniform background with a spike near p=0 where true signal exists. Right: log2FC are centered near 0 with heavier tails for stronger contrasts.</p>
</section>

<section class="card">
  <h2>7. Where is the signal? Cell-type landscape</h2>
  {img(b_heat,"celltype heatmap")}
  <p class="note">Cell type <b>{html.escape(top_ct_overall)}</b> carries the most strong DEGs summed across comparisons. Rows ordered by total signal.</p>
</section>

<h2 style="color:#1f3a5f">8. Per-comparison detail</h2>
{sections}

<footer>Generated by automated EDA · {len(FILES)} files · padj significance = {PADJ_SIG}, strong |log2FC| ≥ {LFC_THR}</footer>
</main></body></html>"""

out = os.path.join(DIR, "differential_expression_EDA.html")
with open(out, "w") as f:
    f.write(HTML)
print("WROTE", out, os.path.getsize(out), "bytes")
