# extract_longest_cds.py

## Description

For each gene in a `*.cds.fa` FASTA file, selects and extracts only the longest CDS isoform. This is useful when a genome annotation produces multiple transcript variants per gene and downstream analyses (e.g., OrthoFinder, phylogenetics) require a single representative sequence per gene.

The script processes all `*.cds.fa` files in the current directory and outputs results to a folder named `最长CDS结果/`.

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
# Run in the directory containing *.cds.fa files
python extract_longest_cds.py
```

No command-line arguments are needed — the script automatically discovers all `*.cds.fa` files in the current working directory.

## Input

FASTA files with headers in the format:

```
>CDS_ID cds ... gene:GENE_ID ...
ATGCCC...
```

The `gene:` tag in the header is used to group CDS isoforms by gene.

## Output (in `最长CDS结果/`)

- `<basename>_cds_info.txt` — tab-separated table of all CDS IDs, gene IDs, and sequence lengths.
- `<basename>_longest_cds.fa` — FASTA file containing only the longest CDS per gene.
