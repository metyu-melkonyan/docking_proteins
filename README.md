# Docking Proteins: Automated Protein-Ligand Docking Analysis

This repository contains Python scripts that automate the analysis of protein-ligand docking results. The workflow focuses on identifying optimal docking solutions, isolating relevant components (such as monomers and ligands), calculating contact frequencies, and preparing the resulting data for further visualization and scientific interpretation.

**Homepage:** [Project Website](https://provartlabundergrads.csb.utoronto.ca/metyus-summative/)

## Features

- **Best Solution Extraction:** Automatically identifies and extracts the best docking poses from docking results.
- **Component Isolation:** Separates key components such as protein monomers and ligand molecules for targeted analysis.
- **Contact Frequency Calculation:** Computes and summarizes contact frequencies between protein and ligand atoms.
- **Data Preparation:** Formats results for downstream visualization tools and interpretation.
- **Automation:** Minimizes manual intervention in the post-docking workflow for rapid and reproducible analysis.

## Getting Started

### Prerequisites

- Python 3.x
- Common scientific libraries: `numpy`, `pandas`, `biopython` (installation commands below)
- Docking output files in a supported format (e.g., PDBQT)

### Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/metyu-melkonyan/docking_proteins.git
    cd docking_proteins
    ```
2. Install dependencies:
    ```bash
    pip install numpy pandas biopython
    ```

### Usage

1. Place your docking result files in the input directory as specified by the script.
2. Run the main analysis script:
    ```bash
    python analyze_docking.py
    ```
3. Output files with extracted solutions, component separations, and summary tables will be generated in the output directory.

*For more detailed usage instructions and options, see comments in the main script or visit the [Project Website](https://provartlabundergrads.csb.utoronto.ca/metyus-summative/).*

## Project Structure

- `analyze_docking.py` – Main automation script for analyzing docking results.
- `utils/` – Helper modules and utility functions.
- `input/` – Place your docking result files here.
- `output/` – Processed results and summary files.

## Citation

If you use this code or workflows in your research, please cite the project or reference the [project website](https://provartlabundergrads.csb.utoronto.ca/metyus-summative/).

## License

This repository currently does not specify a license. Contact the author for permissions regarding usage or distribution.

## Author

[metyu-melkonyan](https://github.com/metyu-melkonyan)
