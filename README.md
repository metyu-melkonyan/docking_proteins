# Protein Docking Analysis Script

## Overview

This script automates the analysis of protein-ligand docking results, focusing on:
- Identifying the best docking poses based on energy scores and residue-ligand interaction frequencies.
- Preparing results for downstream visualization and analysis.
- Normalizing interaction data for comparative studies.

It is designed for computational biology and drug discovery pipelines, offering a robust framework for studying protein-ligand interactions.

---

## Features

- **Best Docking Pose Selection**: Identifies the docking pose with the best energy and interaction scores.
- **Residue-Ligand Interaction Analysis**: Calculates contact frequencies between receptor residues and ligand atoms using distance thresholds.
- **Result Normalization**: Outputs normalized scores for better comparison of docking results.
- **Data Preparation for Visualization**: Converts docking result files into standardized formats compatible with visualization tools like JSmol.

---

## Requirements

- Python 3.8 or higher
- Python libraries:
  - `os`
  - `shutil`
  - `json`
  - `math`
  - `datetime`

---

## Installation

* **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/protein-docking-analysis.git
   cd protein-docking-analysis

# Usage
* Best Docking Pose Selection
Identify the best docking result and prepare it for further analysis:

```python
from docking_analysis import best_result

## Example usage
best_result(
    file_name="file_name.pdb",
    monomer="monomer",
    rec_lig="receptor_ligand",
    receptor="receptor",
    ligand="ligand")
```
* Residue-Ligand Interaction Analysis
Calculate contact frequencies and generate a residue-ligand interaction dictionary:

```python
from docking_analysis import result_dict_generator

## Example usage
result_dict = result_dict_generator(
    threshold=5.0, 
    monomer="monomer", 
    rec_lig="receptor_ligand", 
    receptor="receptor", 
    ligand="ligand")
```
* Visualization-Ready Data Preparation
Normalize interaction scores and prepare files for visualization:

```python 
from docking_analysis import normalize_results

## Example usage
normalize_results(
    monomer_json="monomer_results.json",
    receptor="receptor",
    ligand="ligand",
    final_json="final_normalized_results.json")
```

