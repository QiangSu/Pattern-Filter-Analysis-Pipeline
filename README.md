# scRNA-seq Simulation Pipeline

This repository contains a simulation pipeline to generate Cell Ranger compatible scRNA-seq FASTQ data using Splatter (for counts) and Polyester (for reads), followed by custom processing to add real 10x barcodes, UMIs, and a Poly-T tail.

The pipeline is orchestrated using Snakemake for reproducibility.

## Repository Contents

-   `Snakefile`: Defines the simulation workflow.
-   `config.yaml`: Configurable parameters for the simulation (number of cells, genes, paths, etc.).
-   `environment.yml`: Defines the Conda/Mamba environment with all required software.
-   `scripts/step1_simulate_scRNA_seq_data.R`: R script using Splatter and Polyester.
-   `scripts/step2_add_barcodes_umis_cellranger.py`: Python script to format reads for Cell Ranger.

## Prerequisites

1.  **Git:** To clone this repository.
2.  **Conda or Mamba:** To manage the software environment. Mamba is recommended for speed.
3.  **Reference Transcriptome FASTA:** You will need a transcriptome FASTA file (e.g., from Ensembl or Gencode) for the species you want to simulate. Update the `reference_fasta` path in `config.yaml`.
4.  **10x Barcode Whitelist:** You will need a 10x barcode whitelist file (gzipped). Update the `barcode_whitelist_path` in `config.yaml`.

## Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/your_repo_name.git
    cd your_repo_name
    ```
2.  **Create and activate the Conda/Mamba environment:**
    ```bash
    # If using Mamba:
    mamba env create -f environment.yml
    mamba activate scrna_simulation

    # If using Conda:
    # conda env create -f environment.yml
    # conda activate scrna_simulation
    ```
3.  **Update Configuration:** Edit the `config.yaml` file to set the correct paths to your reference FASTA, barcode whitelist, and customize simulation parameters.

## Running the Pipeline

Execute the Snakemake workflow from the root of the repository:

```bash
# Run locally using the activated environment
snakemake --cores 8 # Use a suitable number of cores

# To run on a cluster, consult Snakemake documentation for cluster execution options
# e.g., using --cluster
```

# FASTQ Base Frequency Logo Generator
This is a Python script designed to visualize the base composition (A, C, G, T, N) at each position within sequencing reads from a FASTQ file. It generates a sequence logo-style plot where the height of each nucleotide character at a given position is proportional to its frequency in the analyzed reads.
**Sequence_logo.py**
## USAGE
```bash
python ./Sequence_logo.py --fastq_file ./T18VN-A30-rep2_score_filtered.fastq.gz --start_pos 1 --end_pos 34 --letter_width 0.8 --output ./T18VN-A30-rep2_score_filtered_1-34posi.png 
```






