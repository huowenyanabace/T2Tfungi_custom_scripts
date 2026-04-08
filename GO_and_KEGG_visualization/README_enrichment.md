# GO & KEGG Enrichment Analysis Visualization Toolkit

Python scripts for generating publication-ready figures from GO and KEGG enrichment analysis results. All figures are output in PDF format at 300 dpi.

## Scripts

### 1. `go_visualization_custom.py` — GO Enrichment Visualization

Generates two figures from GO enrichment results:

- **Bubble plot**: Scatter plot of top enriched GO terms, with bubble size representing gene count, x-axis for fold enrichment, and y-axis for −log10(adjusted p-value). Each term is individually colored.
- **Top terms bar chart**: Horizontal bar chart of top enriched GO terms colored by namespace (BP = red, MF = blue, CC = green).

**Required inputs:**

1. **Enrichment file** (CSV) — must contain columns: `ID`, `Count`, `FoldEnrichment`, `p.adjust`
2. **Description file** (CSV) — must contain column `GO_ID`; optionally `Name` and `Namespace`

**Usage:**

```bash
python3 go_visualization_custom.py -e enrichment.csv -d description.csv -o ./plots -p MyGO
```

| Argument | Flag | Required | Default | Description |
|---|---|---|---|---|
| `--enrichment` | `-e` | Yes | — | GO enrichment results CSV |
| `--description` | `-d` | Yes | — | GO term descriptions CSV |
| `--output_dir` | `-o` | No | `./output` | Output directory |
| `--prefix` | `-p` | No | `GO` | Output filename prefix |

**Outputs:**

- `<prefix>_bubble_plot.pdf`
- `<prefix>_top_terms.pdf`
- `<prefix>_titles.txt` — Figure titles for manuscript use

---

### 2. `kegg_visualization_custom.py` — KEGG Enrichment Visualization

Generates two figures from KEGG enrichment results:

- **Bubble plot**: Scatter plot of top enriched KEGG pathways, with bubble size for gene count, x-axis for fold enrichment, and y-axis for −log10(adjusted p-value).
- **Top pathways bar chart**: Horizontal bar chart colored by fold enrichment intensity (viridis colormap).

**Required input:**

- **Enrichment file** (CSV) — must contain columns: `ID`, `Count`, `FoldEnrichment`, `p.adjust`; optionally `Description`

**Usage:**

```bash
python3 kegg_visualization_custom.py -i enrichment_results.csv -o ./plots -p MyKEGG
```

| Argument | Flag | Required | Default | Description |
|---|---|---|---|---|
| `--input` | `-i` | Yes | — | KEGG enrichment results CSV |
| `--output_dir` | `-o` | No | `./output` | Output directory |
| `--prefix` | `-p` | No | `KEGG` | Output filename prefix |

**Outputs:**

- `<prefix>_bubble_plot.pdf`
- `<prefix>_top_pathways.pdf`
- `<prefix>_titles.txt` — Figure titles for manuscript use

## Requirements

- Python ≥ 3.6
- pandas
- numpy
- matplotlib

```bash
pip3 install pandas numpy matplotlib
```

## Notes

- Both scripts display the top 10 most significant terms/pathways by default (sorted by adjusted p-value).
- Long term or pathway names are automatically truncated in figures.
- The Arial font is used by default. If unavailable, matplotlib falls back to its default font.
- Output directories are created automatically.
