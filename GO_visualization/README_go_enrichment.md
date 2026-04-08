# GO Enrichment Visualization for *Coprinellus radians*

Publication-quality visualization of GO enrichment analysis results for expanded gene families. Available in both Python and R implementations, producing identical figure types.

## Figures Generated

1. **Bubble Plot** — Top 20 enriched GO terms. Bubble size = gene count, color = adjusted p-value, shape = GO category (BP/CC/MF).
2. **GO Hierarchy DAG** — Directed acyclic graph showing hierarchical relationships among enriched GO terms (top 25) and their parental terms. Nodes use numbered labels; a separate legend file maps numbers to GO IDs and term names.

## Requirements

### Python version

- Python ≥ 3.6
- pandas, numpy, matplotlib, networkx

```bash
pip3 install pandas numpy matplotlib networkx
```

### R version

- R ≥ 4.0
- ggplot2, dplyr, tidyr, scales, stringr, igraph, ggraph

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "scales", "stringr", "igraph", "ggraph"))
```

## Input Files

Both scripts require two CSV files:

1. **Enrichment results CSV** — columns must include: `ID`, `Count`, `FoldEnrichment`, `RichFactor`, `GeneRatio`, `p.adjust`
2. **GO annotation CSV** — columns must include: `GO_ID`, `Name`, `Namespace`, `Status`, `Parents`

## Usage

### Python

```bash
python3 plot_go_enrichment.py <enrichment.csv> <go_annotation.csv> [output_prefix]

# Example
python3 plot_go_enrichment.py Crad_GO_enrichment_results.csv Crad_results.csv Crad_GO
```

### R

```bash
Rscript plot_go_enrichment.R <enrichment.csv> <go_annotation.csv> [output_prefix]

# Example
Rscript plot_go_enrichment.R Crad_GO_enrichment_results.csv Crad_results.csv Crad_GO
```

The `output_prefix` argument is optional and defaults to `Crad_GO`.

## Output Files

| File | Description |
|---|---|
| `<prefix>_bubble_plot.pdf / .png` | Bubble plot of top 20 enriched GO terms |
| `<prefix>_GO_DAG.pdf / .png` | GO hierarchy DAG |
| `<prefix>_DAG_node_legend.txt` | Mapping of DAG node numbers to GO IDs and term names |
| `<prefix>_figure_titles.txt` | Figure titles and captions for manuscript use |

## Notes

- All figures are saved at 300 dpi with editable text in PDF (fonttype 42).
- The DAG retains the top 25 enriched terms plus their direct parents and grandparents for hierarchical context.
- Titles are not displayed on figures; they are saved separately in the `_figure_titles.txt` file for flexible use in manuscripts.
