# generate_synteny_links.py

## Description

Generates a synteny link file for **Circos** visualization from JCVI/MCScanX output. The script reads a BED file (gene coordinates) and an anchors file (collinear gene pairs grouped into blocks), then computes the genomic span of each synteny block on both genomes and outputs a six-column link file.

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
# With command-line arguments
python generate_synteny_links.py species.bed species_species.anchors output_links.txt

# Without arguments — uses hardcoded default paths (edit the script)
python generate_synteny_links.py
```

## Arguments (positional)

| Position | Description |
|----------|-------------|
| 1 | BED file with gene coordinates (`chrom  start  end  gene_id`) |
| 2 | Anchors file from JCVI (`gene1  gene2`, blocks separated by `###`) |
| 3 | Output link file path |

## Output

A tab-separated file with six columns per line:

```
chr1  start1  end1  chr2  start2  end2
```

Each line represents one synteny block. This file can be used directly as input for Circos `links` tracks or other genome visualization tools.
