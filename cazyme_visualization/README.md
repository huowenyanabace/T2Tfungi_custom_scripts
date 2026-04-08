# CAZyme Comparative Genomics Visualization Toolkit

A set of Python scripts for visualizing CAZyme (Carbohydrate-Active enZymes) gene distribution across fungal species.

## Scripts

### 1. `plot_cazy_stacked.py` — Stacked Percentage Bar Chart

Generates a stacked bar chart showing the proportional composition of six CAZyme classes (AA, CBM, CE, GH, GT, PL) across species. Total gene counts are displayed above each bar, and percentage labels appear inside segments ≥ 4%.

**Input format** (`TSV`):

| Class | Species_A | Species_B | ... |
|-------|-----------|-----------|-----|
| AA    | 12        | 8         | ... |
| GH    | 45        | 30        | ... |

The first column contains the CAZyme class name; remaining columns are gene counts per organism.

**Usage:**

```bash
# Output to screen
python3 plot_cazy_stacked.py Cazyme_class.tsv

# Save to a directory
python3 plot_cazy_stacked.py Cazyme_class.tsv ./results
```

**Outputs:**

- `cazy_stacked_bar.png` — Stacked percentage bar chart (250 dpi)
- `figure_legend_stacked.txt` — Figure legend with raw and percentage tables

---

### 2. `plot_cazyme_param.py` — Heatmap & Fold-Change Bar Chart

Generates two figures from a family-level CAZyme gene count table:

- **Heatmap**: Row-wise Z-score normalized heatmap of key lignocellulose-degrading CAZyme families (GH, CE, AA, CBM), with raw counts displayed in each cell.
- **Fold-change bar chart**: Horizontal bar chart ranking CAZyme families enriched in the target species group compared to others.

**Input format** (`TSV`):

| Family | Cdom | Crad | Cmic | OtherSp1 | ... |
|--------|------|------|------|----------|-----|
| GH5    | 10   | 8    | 12   | 5        | ... |
| AA9    | 6    | 7    | 5    | 2        | ... |

The first column contains the CAZyme family name; remaining columns are gene counts per organism.

**Target species**: By default, the script treats `Cdom`, `Crad`, and `Cmic` as the target group. To change this, edit the `TARGET_SPECIES` list in the script.

**Usage:**

```bash
python3 plot_cazyme_param.py Cazyme.tsv ./results
```

**Outputs:**

- `heatmap.png` — Z-score heatmap (250 dpi)
- `foldchange_bar.png` — Fold-change horizontal bar chart (250 dpi)
- `figure_legends.txt` — Figure titles and descriptions for both charts

**Fold-change filtering criteria:**

- Target group mean > Others mean × 1.3
- Target group mean ≥ 2

## Requirements

- Python ≥ 3.6
- pandas
- numpy
- matplotlib

Install dependencies:

```bash
pip3 install pandas numpy matplotlib
```

## Notes

- Both scripts attempt to use Noto Sans CJK for Chinese character support. If not available, the system default font is used.
- All output figures are saved at 250 dpi with a white background.
- Output directory is created automatically if it does not exist.
