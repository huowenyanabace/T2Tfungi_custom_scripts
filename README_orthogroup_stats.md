# orthogroup_stats.py

## Description

Generates a comprehensive statistical report from OrthoFinder's `Orthogroups_GeneCount.tsv` output, including:

1. Total gene counts per species across all orthogroups.
2. Species-specific gene counts (genes found in only one species).
3. Single-copy orthologue counts per species.
4. Total and core (present in all species) orthogroup counts.
5. Orthogroups shared exclusively among a user-specified set of species.

## Requirements

- Python 3.6+
- pandas

## Usage

```bash
# Find orthogroups shared only between two species
python orthogroup_stats.py -i Orthogroups_GeneCount.tsv -o stats.txt -s Cdom.pep Crad.pep

# Three species
python orthogroup_stats.py -i Orthogroups_GeneCount.tsv -o stats.txt -s Cdom.pep Crad.pep Cdis.pep
```

## Arguments

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | `Orthogroups_GeneCount.tsv` from OrthoFinder | *required* |
| `-o, --output` | Output report file path | `orthogroup_statistics.txt` |
| `-s, --species` | Two or more species column names to find exclusively shared orthogroups | *required* |

## Output

A formatted text report printed to both the console and the specified output file.
