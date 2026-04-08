# get_longest_pep.py

## Description

Selects the longest protein isoform per gene from `*.pep.fa` FASTA files. This is the protein-level counterpart of `extract_longest_cds.py` and is commonly used to prepare non-redundant proteomes for OrthoFinder, BLAST, or CAZyme annotation.

The script processes all `*.pep.fa` files in the current directory and writes results to `longest_pep_results/`.

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
# Run in the directory containing *.pep.fa files
python get_longest_pep.py
```

## Input

FASTA files with headers containing a `gene:GENE_ID` tag:

```
>PROTEIN_ID cds ... gene:GENE_ID ...
MKTLLF...
```

## Output (in `longest_pep_results/`)

- `<basename>.longest.pep.fa` — FASTA with one longest protein per gene.
- `<basename>.longest.info.txt` — tab-separated table of protein ID, gene ID, and sequence length.
