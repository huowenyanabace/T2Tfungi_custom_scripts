# extract_single_copy_trees.py

## Description

Extracts single-copy orthologue gene trees from OrthoFinder output and prepares them as input for **ASTRAL** species-tree inference. The script:

1. Reads a list of single-copy orthogroup IDs.
2. Extracts corresponding gene trees from OrthoFinder's `Resolved_Gene_Trees.txt`.
3. Simplifies tip labels from full gene names (e.g., `Scere_pep_YNL317W`) to species names (e.g., `Scere`).
4. Filters out any trees that still contain duplicate species (indicating hidden paralogy).
5. Writes the cleaned trees in Newick format, one per line, ready for ASTRAL.

## Requirements

- Python 3.6+
- No external packages required.

## Usage

```bash
# Default file names
python extract_single_copy_trees.py

# Custom file names
python extract_single_copy_trees.py \
  Orthogroups_SingleCopyOrthologues.txt \
  Resolved_Gene_Trees.txt \
  Astral_input.tree
```

## Arguments (positional, optional)

| Position | Description | Default |
|----------|-------------|---------|
| 1 | Single-copy orthogroup ID list | `Orthogroups_SingleCopyOrthologues.txt` |
| 2 | Gene trees file from OrthoFinder | `Resolved_Gene_Trees.txt` |
| 3 | Output file for ASTRAL | `Astral_input.tree` |

## Output

- `Astral_input.tree` — one Newick tree per line with species-level tip labels, directly usable by `java -jar astral.jar -i Astral_input.tree`.
