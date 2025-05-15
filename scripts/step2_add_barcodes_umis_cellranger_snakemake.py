#!/usr/bin/env python3
import gzip
import random
import os
import sys

print("Starting Step 3: Convert FASTA, Add REAL Barcodes/UMIs/PolyT to R1, Create Combined Cell Ranger Structure")

# --- 3.1: Define Paths (Read from Environment Variables) ---
print("Reading paths and parameters from environment variables...")
base_dir = os.environ.get("BASE_DIR")
input_dir = os.environ.get("INPUT_DIR") # This should be the full path to the directory from R script
output_dir = os.environ.get("OUTPUT_DIR") # This should be the full path for Cell Ranger output

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# --- User Defined (Read from Environment Variables) ---
cr_sample_name = os.environ.get("CR_SAMPLE_NAME")
barcode_whitelist_path = os.environ.get("BARCODE_WHITELIST_PATH")

# --- 3.2: Define Simulation Parameters (Read from Environment Variables) ---
n_cells = int(os.environ.get("N_CELLS"))

# Cell Ranger Parameters (check consistency with chosen chemistry/whitelist)
barcode_length = int(os.environ.get("BARCODE_LENGTH"))
umi_length = int(os.environ.get("UMI_LENGTH"))
poly_t_length = int(os.environ.get("POLY_T_LENGTH", "0")) # Default to 0 if not set
cr_read1_length = barcode_length + umi_length + poly_t_length
cr_read2_length = int(os.environ.get("CR_READ2_LENGTH"))

bases = 'ATCG'
placeholder_quality_char = 'I' # Phred score 40 placeholder
poly_t_tail_seq = 'T' * poly_t_length # Pre-calculate poly-T string

print(f"Base Directory: {base_dir}")
print(f"Input Directory (from R script): {input_dir}")
print(f"Output Directory (Cell Ranger): {output_dir}")
print(f"Cell Ranger Sample Name: {cr_sample_name}")
print(f"Barcode Whitelist Path: {barcode_whitelist_path}")
print(f"N Cells: {n_cells}")
print(f"Barcode Length: {barcode_length}")
print(f"UMI Length: {umi_length}")
print(f"Poly-T Length: {poly_t_length}")
print(f"Target CR Read 1 Length: {cr_read1_length}")
print(f"Target CR Read 2 Length: {cr_read2_length}")


# --- Read Barcode Whitelist ---
print(f"Reading barcode whitelist from: {barcode_whitelist_path}")
valid_barcodes = []
try:
    # Use gzip.open for reading compressed file ('rt' for read text mode)
    with gzip.open(barcode_whitelist_path, 'rt') as f:
        for line in f:
            barcode = line.strip()
            if barcode: # Ensure non-empty lines
                # Basic validation (optional but good)
                if len(barcode) == barcode_length and all(c in bases for c in barcode):
                     valid_barcodes.append(barcode)
                else:
                     print(f"Warning: Skipping invalid barcode format in whitelist: {barcode}")

    print(f"Read {len(valid_barcodes)} valid barcodes from whitelist.")
    if not valid_barcodes:
        raise ValueError("Whitelist file is empty or contains no valid barcodes.")
    if len(valid_barcodes) < n_cells:
        raise ValueError(f"Whitelist contains fewer barcodes ({len(valid_barcodes)}) than requested cells ({n_cells}). Cannot proceed.")

except FileNotFoundError:
    print(f"ERROR: Barcode whitelist file not found at {barcode_whitelist_path}")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: Failed to read or process barcode whitelist: {e}")
    sys.exit(1)

# --- Sample unique barcodes from the whitelist for our cells ---
print(f"Sampling {n_cells} unique barcodes from the whitelist...")
try:
    cell_barcodes_assigned = random.sample(valid_barcodes, n_cells)
    print(f"Successfully sampled {len(cell_barcodes_assigned)} unique barcodes.")
except ValueError as e:
     # This should be caught by the check above, but as a safeguard
     print(f"ERROR: Could not sample barcodes. {e}")
     sys.exit(1)

# Map cell input index (0 to n_cells-1) to its assigned barcode
# This determines which barcode goes with reads from sample_1, sample_2 etc.
cell_barcode_map = {i: cell_barcodes_assigned[i] for i in range(n_cells)}

# Save barcode mapping (optional but highly recommended for traceability)
barcode_map_file = os.path.join(output_dir, f"{cr_sample_name}_barcode_mapping.txt")
try:
    with open(barcode_map_file, 'w') as f:
        f.write("InputCellIndex\tAssignedBarcode\n")
        for i in range(n_cells):
             # Map 1-based input cell number (like sample_1) to barcode
             f.write(f"{i+1}\t{cell_barcode_map[i]}\n")
    print(f"Saved barcode mapping to {barcode_map_file}")
except IOError as e:
    print(f"Warning: Could not save barcode map file: {e}")


# --- 3.3: Function to Generate Random UMI ---
def generate_umi(length=umi_length):
    """Generates a random UMI sequence."""
    return ''.join(random.choice(bases) for _ in range(length))

# --- 3.4: Process Input Files and Create Combined Cell Ranger Output ---
print(f"\nProcessing {n_cells} input cell FASTA files...")
print(f"Writing combined output to directory: {output_dir}")

# Define the *single pair* of output filenames using the Cell Ranger sample name
cr_r1_outfile_path = os.path.join(output_dir, f"{cr_sample_name}_S1_R1_001.fastq.gz")
cr_r2_outfile_path = os.path.join(output_dir, f"{cr_sample_name}_S1_R2_001.fastq.gz")

processed_input_files = 0
total_reads_written = 0
error_count_files = 0

# Open the *single* output file pair ONCE before the loop
try:
    # Use gzip.open for writing compressed FASTQ files ('wt' for write text mode)
    with gzip.open(cr_r1_outfile_path, 'wt', compresslevel=6) as cr_r1_out_gz, \
         gzip.open(cr_r2_outfile_path, 'wt', compresslevel=6) as cr_r2_out_gz:

        # Iterate through the input cell files (corresponding to sample_1, sample_2, ...)
        # Note: Polyester outputs sample_1.fasta, sample_2.fasta etc. even if renamed to .fastq
        for i in range(n_cells):
            cell_input_index = i # 0-based index for mapping
            sample_number = i + 1 # 1-based index for filenames
            input_sample_id = f"sample_{sample_number}" # Name from R script output

            # Define input filenames from R script output (they are FASTA content)
            # Polyester R1 contains the sequence that becomes Cell Ranger R2 (cDNA)
            input_r1_polyester_path = os.path.join(input_dir, f"{input_sample_id}_1.fasta") # R script renames to .fasta
            # Polyester R2 is not typically used for GEX simulation structure
            input_r2_polyester_path = os.path.join(input_dir, f"{input_sample_id}_2.fasta") # R script renames to .fasta

            # Get the *assigned valid barcode* for this input cell file
            current_cell_barcode = cell_barcode_map[cell_input_index] # The barcode sampled from the whitelist

            # Check if the primary input file exists
            if not os.path.exists(input_r1_polyester_path):
                print(f"Warning: Input file not found: {input_r1_polyester_path}. Skipping input {input_sample_id}.")
                error_count_files += 1
                continue # Skip to the next input cell file

            # Optional check for the paired file (less critical for GEX simulation)
            if not os.path.exists(input_r2_polyester_path):
                 print(f"Warning: Paired input file not found: {input_r2_polyester_path}. Proceeding using only {input_r1_polyester_path}.")

            if sample_number % 500 == 0: # Print progress
                 print(f"  Processing input file {sample_number}/{n_cells} ({input_sample_id})...")

            # Process the reads from this input cell file (FASTA format - 2 lines per record)
            try:
                # Open the *individual* input Polyester R1 FASTA file for this cell
                with open(input_r1_polyester_path, 'r') as poly_r1_in:
                    read_counter_this_file = 0
                    # Process reads 2 lines at a time (one FASTA record)
                    while True:
                        poly_r1_l1_header = poly_r1_in.readline().strip() # Header line starting with >
                        if not poly_r1_l1_header: break # End of this input file
                        poly_r1_l2_seq = poly_r1_in.readline().strip() # Sequence line

                        # --- Create Cell Ranger R1 Record (VALID Barcode + UMI + PolyT) ---
                        umi = generate_umi(umi_length)
                        # Concatenate Barcode, UMI, and Poly-T tail
                        cr_r1_seq = current_cell_barcode + umi + poly_t_tail_seq
                        # Use placeholder quality for R1 (length matches cr_r1_seq)
                        cr_r1_qual = placeholder_quality_char * len(cr_r1_seq)

                        # --- Create Cell Ranger R2 Record (cDNA) ---
                        cr_r2_seq = poly_r1_l2_seq
                        # Ensure R2 sequence length matches target (trim if necessary)
                        if len(cr_r2_seq) > cr_read2_length:
                            cr_r2_seq = cr_r2_seq[:cr_read2_length]
                        # Use placeholder quality for R2 (length matches cr_r2_seq)
                        cr_r2_qual = placeholder_quality_char * len(cr_r2_seq)


                        # --- Construct FASTQ Headers ---
                        # Make header unique across *all* reads in the output files
                        # Use original header info to make read_name unique and trace back
                        # Example: @<Polyester_ID> 1:N:0:BARCODE
                        # Let's extract the ID from the Polyester header (after '>')
                        poly_read_id = poly_r1_l1_header.lstrip('>').split(' ')[0] # Get ID before first space
                        # Use sample_number and original read ID to form unique ID
                        cr_header_base = f"@{sample_number}_{poly_read_id}"
                        cr_r1_header = f"{cr_header_base} 1:N:0:{current_cell_barcode}" # Read 1, No filter, Control 0, Barcode
                        cr_r2_header = f"{cr_header_base} 2:N:0:{current_cell_barcode}" # Read 2, No filter, Control 0, Barcode


                        # --- Write to the SINGLE combined output files ---
                        cr_r1_out_gz.write(f"{cr_r1_header}\n{cr_r1_seq}\n+\n{cr_r1_qual}\n")
                        cr_r2_out_gz.write(f"{cr_r2_header}\n{cr_r2_seq}\n+\n{cr_r2_qual}\n")

                        read_counter_this_file += 1
                        total_reads_written += 1

                processed_input_files += 1

            except Exception as e:
                print(f"Error processing reads from input file {input_sample_id}: {e}")
                # Optionally add to a separate error counter for read processing errors
                continue # Move to the next input file

except Exception as e:
    print(f"FATAL ERROR: Could not open or write to output files in {output_dir}. Check permissions and path.")
    print(f"Error details: {e}")
    # Clean up potentially partially written files if error occurs during opening/writing
    if os.path.exists(cr_r1_outfile_path): os.remove(cr_r1_outfile_path)
    if os.path.exists(cr_r2_outfile_path): os.remove(cr_r2_outfile_path)
    sys.exit(1)


# --- 3.5: Final Summary ---
print("\n--- Processing Summary ---")
print(f"Total input cell FASTA files targeted: {n_cells}")
print(f"Successfully processed input files: {processed_input_files}")
print(f"Input files skipped due to errors/missing: {error_count_files}")
print(f"Total reads written to combined output files: {total_reads_written}")
if error_count_files > 0:
    print("Warning: Some input files were skipped. Check logs above.")

if processed_input_files > 0:
    print(f"Step 3: Created combined Cell Ranger compatible FASTQ files using whitelist and added Poly-T tail.")
    print(f"  Output R1: {cr_r1_outfile_path} (Length: {cr_read1_length} bp = {barcode_length}B + {umi_length}U + {poly_t_length}T)")
    print(f"  Output R2: {cr_r2_outfile_path} (Length: max {cr_read2_length} bp)")
    print(f"  Barcode whitelist used: {barcode_whitelist_path}")
    print(f"  Barcode-to-Input mapping saved: {barcode_map_file}")
else:
    print("\nERROR: No input files were successfully processed. Check input directory and naming.")
    sys.exit(1)
