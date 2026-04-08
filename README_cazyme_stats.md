# cazyme_stats.py

## Description

Summarizes CAZyme (Carbohydrate-Active Enzyme) annotation results produced by **dbCAN2 / run_dbcan**. The script scans all `*_result` directories in a working directory, parses each `overview.txt` file, and generates two tab-separated statistics tables:

1. **Major-class table** — gene counts for each of the six CAZyme classes (GH, GT, PL, CE, AA, CBM) per species.
2. **Subfamily table** — gene counts for every individual CAZyme family (e.g., GH5, AA9) per species.

Only annotations supported by at least *N* tools (default 2 of 3: HMMER, DIAMOND, Hotpep) are retained, following standard dbCAN2 best-practice filtering.

## Requirements

- Python 3.6+
- No external packages required (uses only the standard library).

## Usage

```bash
# Basic usage — scan current directory
python cazyme_stats.py -d . -o cazyme_stats

# Require at least 2 supporting tools (default)
python cazyme_stats.py -d /path/to/work -o cazyme_stats -m 2

# Custom result-folder pattern
python cazyme_stats.py -d . -o cazyme_stats -p "*_dbcan"
```

## Arguments

| Flag | Description | Default |
|------|-------------|---------|
| `-d, --dir` | Working directory containing result folders | `.` |
| `-o, --output` | Output file prefix | `cazyme_stats` |
| `-m, --min-tools` | Minimum number of tools that must agree | `2` |
| `-p, --pattern` | Glob pattern for result folders | `*_result` |

## Input

Each result folder (e.g., `Cdis_result/`) must contain an `overview.txt` file produced by `run_dbcan`.

## Output

- `<prefix>_major_classes.tsv` — per-species counts for GH, GT, PL, CE, AA, CBM, and totals.
- `<prefix>_subfamilies.tsv` — per-species counts for every CAZyme family detected.
- Console summary table.
