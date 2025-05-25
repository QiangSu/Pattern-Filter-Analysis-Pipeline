# Scanpy Single-Cell Analysis Project

This directory contains a dedicated pipeline for single-cell RNA sequencing (scRNA-seq) data analysis using the Scanpy framework. It integrates custom cell typing models and reference gene markers to facilitate in-depth biological insights, particularly focusing on comparative analysis and potentially cell type ratio dynamics (as suggested by the notebook name).

## Table of Contents

-   [Overview](#overview)
-   [Features](#features)
-   [Project Structure](#project-structure)
-   [Setup and Installation](#setup-and-installation)
    -   [Prerequisites](#prerequisites)
    -   [Environment Setup](#environment-setup)
    -   [Data and Model Preparation](#data-and-model-preparation)
-   [Usage](#usage)
-   [Results](#results)
-   [Notes and Considerations](#notes-and-considerations)
-   [Contact](#contact)
-   [License](#license)

## Overview

This project provides a robust workflow for analyzing scRNA-seq datasets. It leverages the power of Scanpy for processing, visualization, and interpretation of single-cell data. The primary analysis is conducted within the `scanpy_analysis_ratio_git.ipynb` Jupyter notebook, guiding users through data loading, quality control, normalization, dimensionality reduction, clustering, and potentially cell type annotation using an integrated CellTypist model and custom marker references.

## Features

-   **Comprehensive scRNA-seq Analysis:** Full pipeline from raw data to processed AnnData objects and insightful visualizations.
-   **Cell Type Annotation:** Utilizes a pre-trained `Mouse_Whole_Brain.pkl` CellTypist model for automated cell type identification.
-   **Reference Marker Integration:** Incorporates `cellmarker_mouse_brain_seq_exp_review.csv` for validating cell types or performing targeted differential expression analysis.
-   **Modular Design:** Data, models, references, and results are logically separated for clarity and reproducibility.
-   **Jupyter Notebook Workflow:** Interactive and step-by-step analysis with explanations.

## Project Structure

scanpy_analysis_project/
├── scanpy_analysis_ratio_git.ipynb    # Main Jupyter Notebook for the analysis workflow
├── data/
│   └── filtered_feature_bc_matrix_WT_filtereddata_hd1/ # Input raw scRNA-seq data (e.g., 10x Genomics output)
├── models/
│   └── Mouse_Whole_Brain.pkl          # Pre-trained CellTypist model for cell type annotation
├── references/
│   └── cellmarker_mouse_brain_seq_exp_review.csv # Curated list of mouse brain cell markers
├── results/                           # Directory to store analysis outputs (plots, tables, AnnData objects)
├── .gitignore                         # Local Git ignore rules (e.g., for ‘results/’ and intermediate files)
└── requirements.txt                   # Python dependencies specific to this Scanpy analysis


## Setup and Installation

This `scanpy_analysis_project` is designed to be a subdirectory within the larger `Pattern-Filter-Analysis-Pipeline` repository.

### Prerequisites

-   Git
-   Conda (Miniconda or Anaconda recommended for managing Python environments)

### Environment Setup

1.  **Clone the main repository:**
    If you haven't already, clone the parent `Pattern-Filter-Analysis-Pipeline` repository:
    ```bash
    git clone https://github.com/QiangSu/Pattern-Filter-Analysis-Pipeline.git
    cd Pattern-Filter-Analysis-Pipeline
    ```

2.  **Navigate to the project directory:**
    ```bash
    cd scanpy_analysis_project
    ```

3.  **Create and activate a Conda environment:**
    It's highly recommended to create a dedicated Conda environment to manage dependencies for this project.
    ```bash
    conda create -n scanpy_env python=3.9  # Or your preferred Python version (e.g., 3.8, 3.10)
    conda activate scanpy_env
    ```

4.  **Install required Python packages:**
    The `requirements.txt` file lists all necessary libraries.
    ```bash
    pip install -r requirements.txt
    # If there's a consolidated requirements.txt at the root of Pattern-Filter-Analysis-Pipeline,
    # you might install from there instead: pip install -r ../requirements.txt
    ```
    *Note: If you encounter issues with `scanpy` or `anndata` installation, sometimes installing from `conda-forge` channels helps:*
    `conda install -c conda-forge scanpy jupyterlab`

### Data and Model Preparation

The `data/`, `models/`, and `references/` directories are pre-populated with example files as part of this repository.

-   **`data/filtered_feature_bc_matrix_WT_filtereddata_hd1/`**: This directory is expected to contain raw 10x Genomics data in the standard Cell Ranger output format (e.g., `matrix.mtx`, `features.tsv`, `barcodes.tsv`). Ensure your actual experimental data for analysis is placed here following this structure.
-   **`models/Mouse_Whole_Brain.pkl`**: This is a pre-trained CellTypist model. If you wish to use a different model, replace this file with your desired `.pkl` model.
-   **`references/cellmarker_mouse_brain_seq_exp_review.csv`**: This CSV file contains reference cell markers. You can modify or replace this file with your own marker gene lists as needed.

## Usage

The primary analysis workflow is encapsulated in the `scanpy_analysis_ratio_git.ipynb` Jupyter notebook.

1.  **Activate your Conda environment:**
    ```bash
    conda activate scanpy_env
    ```

2.  **Launch Jupyter Lab (or Jupyter Notebook):**
    Ensure you are in the `scanpy_analysis_project/` directory in your terminal before launching.
    ```bash
    jupyter lab
    ```
    This will open Jupyter Lab in your web browser.

3.  **Open the Notebook:**
    Navigate to and open `scanpy_analysis_ratio_git.ipynb`.

4.  **Execute Cells:**
    Run the cells sequentially within the notebook. The notebook is structured to guide you through:
    -   Loading the scRNA-seq data.
    -   Initial quality control and filtering.
    -   Normalization and scaling.
    -   Dimensionality reduction (PCA, UMAP).
    -   Clustering.
    -   Cell type annotation using CellTypist and marker gene validation.
    -   Any specific "ratio" analysis or comparative studies defined within the notebook.

## Results

The `results/` directory will be populated with various outputs generated during the analysis, including:

-   High-resolution plots (e.g., UMAPs, gene expression heatmaps, violin plots)
-   Tables of differentially expressed genes
-   Annotated AnnData objects in `.h5ad` format for downstream analysis
-   Any custom output specific to the "ratio" analysis performed.

This directory is typically ignored by Git (as per `.gitignore`) to prevent committing large data products.

## Notes and Considerations

-   **Computational Resources:** scRNA-seq analysis can be memory-intensive. Ensure your system has sufficient RAM (e.g., 16GB+ for typical datasets; 32GB+ for larger datasets).
-   **Data Format:** The notebook expects 10x Genomics raw data output in the `data/` directory. If your data is in a different format (e.g., a pre-processed AnnData object), you may need to adjust the data loading steps in the notebook.
-   **Customization:** Feel free to modify the `Mouse_Whole_Brain.pkl` model or `cellmarker_mouse_brain_seq_exp_review.csv` reference file to suit your specific experimental context or species.
-   **Jupyter Kernel:** Ensure the Jupyter notebook is running with the `scanpy_env` kernel. You can select this from the "Kernel" menu in Jupyter Lab.

## Contact

For questions or issues related to this specific Scanpy analysis project, please contact:

Qiang Su
[qiang.su@example.com](mailto:qiang.su@example.com) (Please replace with your actual contact email)

## License

This project is part of the `Pattern-Filter-Analysis-Pipeline` and is subject to its overarching license (if specified, typically MIT or Apache 2.0).

