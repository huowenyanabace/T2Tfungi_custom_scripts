# extract_genes_from_orthogroups.py

## Description

Extracts gene IDs for a specific species from selected orthogroups in an OrthoFinder `Orthogroups.tsv` file. Given a list of orthogroup IDs (e.g., significantly expanded gene families from CAFE), the script pulls out all gene IDs belonging to a target species column and writes them to a plain-text file (one gene per line). This gene list can then be used as input for downstream enrichment analyses.

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
# Basic usage
python extract_genes_from_orthogroups.py \
  -i expanded_OG_IDs.txt \
  -o Orthogroups.tsv \
  -c Cdom.pep \
  -out expanded_genes.txt

# Show all available column names
python extract_genes_from_orthogroups.py \
  -o Orthogroups.tsv --show-columns
```

## Arguments

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --id-file` | Text file with one orthogroup ID per line (e.g., `OG0000123`) | *required* |
| `-o, --orthogroups` | OrthoFinder `Orthogroups.tsv` file | *required* |
| `-c, --column` | Column name of the target species (e.g., `Cdom.pep`) | *required* |
| `-out, --output` | Output file path | `extracted_genes.txt` |
| `--show-columns` | Print all column names and exit | — |

## Output

A plain-text file with one gene ID per line (deduplicated), suitable for clusterProfiler or similar tools.
