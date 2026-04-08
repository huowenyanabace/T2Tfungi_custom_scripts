# eggnog_enrichment.R

## Description

Performs GO and KEGG enrichment analysis on a list of expanded gene families using eggNOG-mapper annotation results and **clusterProfiler**. The script:

1. Parses an eggNOG-mapper `.annotations` file (auto-fixes column-count mismatches).
2. Builds TERM2GENE mappings for GO and KEGG.
3. Reads a target gene list from `expanded_genes.txt`.
4. Runs over-representation analysis (`enricher`) for both GO and KEGG.
5. Generates bar plots, dot plots, cnet plots, and heat plots as PDF.

## Requirements

- R 4.0+
- R packages: `clusterProfiler`, `enrichplot`, `tidyr`, `dplyr`, `ggplot2`

## Usage

Before running, edit the script to set:

- `eggnog_file` — path to your `.emapper.annotations` file.
- `expanded_genes.txt` — one gene ID per line (or define genes in-script).

```bash
Rscript eggnog_enrichment.R
```

## Input

- eggNOG-mapper annotation file (`.emapper.annotations`)
- `expanded_genes.txt` — gene list of interest (one ID per line)

## Output

- `GO_enrichment_results.csv`, `KEGG_enrichment_results.csv`
- `GO_barplot.pdf`, `GO_dotplot.pdf`, `GO_cnetplot.pdf`, `GO_heatplot.pdf`
- `KEGG_barplot.pdf`, `KEGG_dotplot.pdf`, `KEGG_cnetplot.pdf`
- `enrichment_analysis_summary.csv`
- `eggnog_annotation_data.RData`
