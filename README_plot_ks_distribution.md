# plot_ks_distribution.R

## Description

Plots synonymous substitution rate (Ks) distribution curves for detecting whole-genome duplication (WGD) events. The script reads a tab-separated file containing Ks values for paralogous and/or orthologous gene pairs, bins them, and draws frequency distribution curves colored by species or comparison type.

## Requirements

- R 4.0+
- R packages: `ggplot2`

## Usage

Edit the following variables in the script before running:

- `input_file` — path to the Ks data file (default: `all.ks`).
- `output_file` — output PDF filename (default: `Ks_distribution.pdf`).
- `colors` — named color vector matching your species/comparison labels.

```bash
Rscript plot_ks_distribution.R
```

## Input

A tab-separated file with three columns (header required):

```
pair        Ks      Species
gene1-gene2 0.523   Cdis
gene3-gene4 1.204   Cdis_Cdom
```

- **pair** — gene pair identifier
- **Ks** — synonymous substitution rate
- **Species** — species label or pairwise comparison label

## Output

- `Ks_distribution.pdf` — line plot of Ks frequency distributions.
- `Ks_distribution.png` — same plot in PNG format.
- Console summary statistics (mean, median, range) per species/comparison.

## Parameters

| Variable | Description | Default |
|----------|-------------|---------|
| `ks_max` | Maximum Ks value to include | `5` |
| `bin_width` | Histogram bin width | `0.1` |
