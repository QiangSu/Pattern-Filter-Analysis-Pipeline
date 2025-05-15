# Snakefile
# Workflow for simulating scRNA-seq data compatible with Cell Ranger

# Load configuration
configfile: "config.yaml"

# Define base directories and full paths using config
BASE_DIR = config['base_dir']
OUTPUT_FASTQ_DIR = os.path.join(BASE_DIR, config['output_fastq_dir_rel'])
CELLRANGER_OUTPUT_DIR = os.path.join(BASE_DIR, config['cellranger_output_dir_rel'])

# Define full paths for specific output files
GROUND_TRUTH_MATRIX = os.path.join(BASE_DIR, config['ground_truth_matrix_rel'])
GROUND_TRUTH_CELLTYPES = os.path.join(BASE_DIR, config['ground_truth_celltypes_rel'])
CR_R1_OUT = os.path.join(CELLRANGER_OUTPUT_DIR, config['cr_r1_output'])
CR_R2_OUT = os.path.join(CELLRANGER_OUTPUT_DIR, config['cr_r2_output'])
BARCODE_MAP_FILE = os.path.join(CELLRANGER_OUTPUT_DIR, config['barcode_map_file'])


# --- Rule all: Defines the final target files ---
# Running 'snakemake' without a specific target will build these files
rule all:
    input:
        CR_R1_OUT,
        CR_R2_OUT,
        GROUND_TRUTH_MATRIX,
        GROUND_TRUTH_CELLTYPES,
        BARCODE_MAP_FILE

# --- Rule 1: Simulate Counts and Initial FASTA files using R script ---
# This rule produces the ground truth files and the directory containing per-cell FASTAs
rule simulate_counts_and_initial_fasta:
    input:
        script = config['scripts']['simulate_r'],
        reference_fasta = config['reference_fasta']
    output:
        # Depend on key outputs from the R script
        directory(OUTPUT_FASTQ_DIR), # Declare the output directory
        GROUND_TRUTH_MATRIX,
        GROUND_TRUTH_CELLTYPES
    params:
        # Pass parameters from config to the script via environment variables
        base_dir = BASE_DIR,
        output_fastq_dir_rel = config['output_fastq_dir_rel'],
        ground_truth_matrix_rel = config['ground_truth_matrix_rel'],
        ground_truth_celltypes_rel = config['ground_truth_celltypes_rel'],
        reference_fasta = config['reference_fasta'],
        n_cells = config['n_cells'],
        n_genes = config['n_genes'],
        n_groups = config['n_groups'],
        cr_read2_length = config['cr_read2_length'] # R script uses this for readlen
    shell:
        """
        mkdir -p {params.base_dir}
        mkdir -p {OUTPUT_FASTQ_DIR}

        BASE_DIR={params.base_dir} \\
        OUTPUT_FASTQ_DIR_REL={params.output_fastq_dir_rel} \\
        GROUND_TRUTH_MATRIX_REL={params.ground_truth_matrix_rel} \\
        GROUND_TRUTH_CELLTYPES_REL={params.ground_truth_celltypes_rel} \\
        REFERENCE_FASTA={params.reference_fasta} \\
        N_CELLS={params.n_cells} \\
        N_GENES={params.n_genes} \\
        N_GROUPS={params.n_groups} \\
        CR_READ2_LENGTH={params.cr_read2_length} \\
        Rscript {input.script}
        """

# --- Rule 2: Convert FASTA to Cell Ranger FASTQ using Python script ---
# This rule takes the per-cell FASTAs and barcode whitelist,
# adds barcodes/UMIs/PolyT, and combines into the final gzipped FASTQ pair.
rule create_cellranger_fastq:
    input:
        script = config['scripts']['process_python'],
        # Depend on the output directory from the R script rule
        polyester_fastq_dir = OUTPUT_FASTQ_DIR,
        barcode_whitelist = config['barcode_whitelist_path']
    output:
        # The final gzipped R1 and R2 files for Cell Ranger
        CR_R1_OUT,
        CR_R2_OUT,
        BARCODE_MAP_FILE # Also output the barcode mapping file
    params:
        # Pass parameters from config to the script via environment variables
        base_dir = BASE_DIR,
        input_dir = OUTPUT_FASTQ_DIR, # Pass the full path
        output_dir = CELLRANGER_OUTPUT_DIR, # Pass the full path
        cr_sample_name = config['cr_sample_name'],
        barcode_whitelist_path = config['barcode_whitelist_path'],
        n_cells = config['n_cells'],
        barcode_length = config['barcode_length'],
        umi_length = config['umi_length'],
        poly_t_length = config['poly_t_length'],
        cr_read2_length = config['cr_read2_length']
    shell:
        """
        mkdir -p {CELLRANGER_OUTPUT_DIR}

        BASE_DIR={params.base_dir} \\
        INPUT_DIR={params.input_dir} \\
        OUTPUT_DIR={params.output_dir} \\
        CR_SAMPLE_NAME={params.cr_sample_name} \\
        BARCODE_WHITELIST_PATH={params.barcode_whitelist_path} \\
        N_CELLS={params.n_cells} \\
        BARCODE_LENGTH={params.barcode_length} \\
        UMI_LENGTH={params.umi_length} \\
        POLY_T_LENGTH={params.poly_t_length} \\
        CR_READ2_LENGTH={params.cr_read2_length} \\
        python {input.script}
        """

# Define how to create the directories needed before rules run
# Snakemake will automatically create output directories, but explicitly defining
# the base_dir can be helpful if it's not an output of any rule.
# rule create_base_dir:
#    output: directory(BASE_DIR)
#    shell: "mkdir -p {output}"

# The rules implicitly create directories specified in their outputs.
# The shell commands in the rules also include `mkdir -p` as a failsafe.


