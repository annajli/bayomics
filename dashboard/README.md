# SoundLife Capstone — EDA Dashboard

Interactive Streamlit dashboard that unifies the exploratory data analysis for
the SoundLife healthy-cohort capstone. It pulls together the four EDA outputs
in `../eda/` plus the live DESeq2 result tables.

## Pages

| Page | Source | Type |
|---|---|---|
| **Overview** | summary metrics + contrast key | KPIs |
| **Differential Expression** | 4 DESeq2 result CSVs | 🔵 Interactive (volcano, filters, gene search, cell-type ranking) |
| **Clinical Labs & Metadata** | `../eda/clinical_labs_metadata.ipynb` | Rendered notebook (plots + narrative) |
| **Flu Serology & HAI** | `../eda/flu_serology_hai.ipynb` | Rendered notebook (plots + narrative) |
| **Pseudobulk Report** | `../eda/EDA_pseudobulk_report.html` | Embedded report |
| **Plasma & Whole Blood** | `../eda/plasma_wholeblood.html` | Embedded report |

## Setup

```bash
cd dashboard
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

Opens at <http://localhost:8501>.

## Layout

```
dashboard/
├── app.py            # the dashboard (multi-page via st.navigation)
├── requirements.txt
├── README.md
└── assets/
    ├── clinical_labs_metadata.html   # pre-rendered notebook
    └── flu_serology_hai.html         # pre-rendered notebook
```

## Notes

- **Live data:** the **Differential Expression** page reads the raw DESeq2 CSVs
  from `../../UVA_Data_Capstone/Datasets/Differential_Expression/`. If those
  files are missing, the page automatically falls back to the static
  `../eda/differential_expression_EDA.html` report.
- **Rendered notebooks:** the clinical-labs and flu-serology pages embed HTML
  exported from the notebooks (using their already-saved plot outputs; code
  cells hidden). They are **not** re-executed, so they don't need the raw
  serology/labs CSVs.
- **Regenerating the notebook HTML** (only needed if you re-run the notebooks):

  ```bash
  jupyter nbconvert --to html --template classic --no-input \
    ../eda/clinical_labs_metadata.ipynb --output-dir assets
  jupyter nbconvert --to html --template classic --no-input \
    ../eda/flu_serology_hai.ipynb --output-dir assets
  ```

- Thresholds used throughout: significance `padj < 0.05`; "strong" hit
  additionally requires `|log2FC| ≥ 1`. Adjust `PADJ_SIG` / `LFC_THR` at the top
  of `app.py`.
