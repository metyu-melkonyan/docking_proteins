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

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/protein-docking-analysis.git
   cd protein-docking-analysis
