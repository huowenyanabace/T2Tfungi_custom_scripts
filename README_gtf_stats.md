# gtf_stats.py

## Description

Computes basic structural statistics from a GTF gene annotation file:

- Total number of genes
- Average gene length (bp)
- Average number of exons per gene
- Average intron length (bp)
- Total exon and intron counts

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
python gtf_stats.py annotation.gtf
```

## Arguments

| Position | Description |
|----------|-------------|
| 1 | Path to the GTF file | 

## Output

Statistics are printed to standard output. Example:

```
Total genes:           10324
Average gene length:   2156.43 bp
Exons per gene:        5.50
Average intron length: 112.87 bp
```
