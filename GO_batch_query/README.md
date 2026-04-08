# GO Term Batch Query Tool User Guide

## Features
- Batch reading of GO IDs from a txt file
- Two output formats: text and CSV
- Automatic handling of blank lines and comment lines
- Detailed query results, including definitions, parent nodes, child nodes, etc.

## Installation

### 1. Install Dependencies
```bash
# Update system
sudo apt update

# Install Python3 and pip
sudo apt install python3 python3-pip -y

# Install goatools
pip3 install goatools
```

### 2. Download the GO Database
```bash
# Download go-basic.obo (recommended, ~40MB)
wget http://purl.obolibrary.org/obo/go/go-basic.obo

# Or download the full go.obo (larger, includes more relationships)
# wget http://purl.obolibrary.org/obo/go.obo
```

## Usage

### 1. Prepare the GO ID File
Create a txt file (e.g., `go_ids.txt`) with one GO ID per line:
```
GO:0008150
GO:0005515
GO:0005575
```

Comments with `#` are supported:
```
# This is a comment line
GO:0008150  # Biological process
GO:0005515  # Protein binding
```

### 2. Run a Query

**Basic usage (output to screen):**
```bash
python3 go_batch_query.py go_ids.txt
```

**Save as a text file:**
```bash
python3 go_batch_query.py go_ids.txt -o results.txt
```

**Save as a CSV file:**
```bash
python3 go_batch_query.py go_ids.txt -o results.csv -f csv
```

**Use a custom OBO file:**
```bash
python3 go_batch_query.py go_ids.txt --obo /path/to/go.obo -o results.txt
```

### 3. View Help
```bash
python3 go_batch_query.py -h
```

## Output Format

### Text Output
```
======================================================================
[1/5] GO ID: GO:0008150
======================================================================
Name: biological_process
Namespace: biological_process
Definition: A biological process is the execution of a molecular-level...
Depth: 0
Parents: None (root node)
Children Count: 31
======================================================================
```

### CSV Output
Contains the following columns:
- GO_ID: GO identifier
- Name: GO term name
- Namespace: Namespace (biological_process / molecular_function / cellular_component)
- Definition: Definition
- Depth: Hierarchy depth
- Parents: Parent node information
- Children_Count: Number of child nodes
- Status: Query status (Found / Not Found)

## Example Files
- `go_batch_query.py`: Main script
- `go_ids_example.txt`: Example GO ID file

## FAQ

**Q: "go-basic.obo file not found" error?**
A: Make sure the `go-basic.obo` file is in the same directory where you run the script, or specify the file path using the `--obo` argument.

**Q: Some GO IDs return "NOT FOUND"?**
A: Check that the GO ID format is correct (should be `GO:xxxxxxx`). The GO ID may also have been deprecated.

**Q: How do I update the GO database?**
A: Simply re-download the latest `go-basic.obo` file:
```bash
wget http://purl.obolibrary.org/obo/go/go-basic.obo
```

## Performance Notes
- The initial loading of the GO database takes 10–30 seconds.
- After that, each GO ID query completes almost instantly.
- Batch querying is recommended to avoid reloading the database repeatedly.

## Troubleshooting
If you encounter issues, please check:
1. Python version >= 3.6
2. goatools is properly installed
3. The GO database file exists and is readable
