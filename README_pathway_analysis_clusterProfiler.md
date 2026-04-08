# pathway_analysis_clusterProfiler.R

## Description

Performs a genome-wide functional annotation survey from an eggNOG-mapper output file. Unlike `eggnog_enrichment.R` (which tests a specific gene list), this script characterizes the **entire genome's** annotation landscape:

1. Counts and visualizes GO term distribution (top 30 bar chart).
2. Counts and visualizes KEGG KO and pathway distributions.
3. Counts and visualizes COG functional category distribution (bar chart + pie chart).
4. Generates an annotation coverage summary (percentage of genes with GO/KEGG/COG).
5. Exports TERM2GENE mapping files for downstream enrichment analyses.

## Requirements

- R 4.0+
- R packages: `tidyverse`, `clusterProfiler`, `enrichplot`, `ggplot2`, `scales`

## Usage

Edit the `emapper_file` variable in the script to point to your annotation file, then run:

```bash
Rscript pathway_analysis_clusterProfiler.R
```

## Input

- eggNOG-mapper annotation file (`.emapper.annotations`)

## Output

- `GO_term_statistics.txt`, `GO_term_distribution.pdf`
- `KEGG_KO_statistics.txt`, `KEGG_pathway_statistics.txt`, `KEGG_pathway_distribution.pdf`
- `COG_category_statistics.txt`, `COG_category_distribution.pdf`, `COG_category_pie_chart.pdf`
- `annotation_summary.txt`, `annotation_coverage.pdf`
- `GO_TERM2GENE.txt`, `KEGG_KO_TERM2GENE.txt`, `KEGG_PATHWAY_TERM2GENE.txt`, `COG_TERM2GENE.txt`
