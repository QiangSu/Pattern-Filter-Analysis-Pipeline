#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import gzip

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a base frequency sequence logo from a FASTQ file.')
parser.add_argument('--fastq_file', type=str, required=True, help='Path to the input FASTQ file (.fastq or .fastq.gz)')
parser.add_argument('--start_pos', type=int, default=1, help='Starting position for plotting (1-based, default: 1)')
parser.add_argument('--end_pos', type=int, default=150, help='Ending position for plotting (1-based, default: 150)')
parser.add_argument('--output', type=str, help='Output file path to save the plot (e.g., plot.png or plot.pdf)')
parser.add_argument('--letter_width', type=float, default=0.8, help='Width of each letter in the plot (default: 0.8)')
args = parser.parse_args()

# Use the provided arguments
fastq_file = args.fastq_file
# Convert start_pos and end_pos to 0-based indices IMMEDIATELY
start_pos = args.start_pos - 1
end_pos = args.end_pos - 1
output_file = args.output
letter_width = args.letter_width


# Dictionary to hold counts. Key = position index, value = dict of base counts
positions = {}
total_reads = 0  # Variable to count total number of reads

# Function to open FASTQ file, handling both .fastq and .fastq.gz
def open_fastq(filepath):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt")  # "rt" mode for text reading
    else:
        return open(filepath, "r")

# Process the FASTQ file and count reads
with open_fastq(fastq_file) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        total_reads += 1  # Increment read count
        seq_str = str(record.seq)
        for i, base in enumerate(seq_str):
            if i not in positions:
                positions[i] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
            positions[i][base] += 1

# Convert counts to frequencies
freqs_data = {}
for i in positions:
    total_count = sum(positions[i].values())
    freqs_data[i] = {}
    for base in ["A", "C", "G", "T", "N"]:
        freqs_data[i][base] = positions[i][base] / total_count if total_count > 0 else 0

# Prepare data for plotting
plot_positions = sorted(positions.keys())
num_plot_positions = len(plot_positions)

# Plot frequencies as stretched alphabet characters
plt.figure(figsize=(40, 15))  # Custom width and height in inches
ax = plt.gca()  # Get current axes

# Use the 0-based start_pos and end_pos for filtering
plot_positions_range = [p for p in plot_positions if start_pos <= p <= end_pos]
num_plot_positions_range = len(plot_positions_range)

if not plot_positions_range:
    print(f"No positions to plot within the range {start_pos + 1}-{end_pos + 1}. Adjust start_pos and end_pos.") # Display 1-based range in error message
    sys.exit()

x_coords = np.arange(num_plot_positions_range)
# Generate 1-based labels AFTER filtering
plot_labels = [pos + 1 for pos in plot_positions_range]
plt.xticks(x_coords, plot_labels)

# Make tick labels bold and set size
plt.tick_params(axis='x', labelsize=25, width=2, length=6, which='major')  # Adjust width and length as needed
plt.tick_params(axis='y', labelsize=25, width=2, length=6, which='major')
for tick in ax.get_xticklabels():
    tick.set_fontweight("bold")
for tick in ax.get_yticklabels():
    tick.set_fontweight("bold")


# Set x and y labels with bold font
plt.xlabel('Position in Read', fontsize=30, fontweight='bold')
plt.ylabel('Frequency', fontsize=30, fontweight='bold')

# Update title to include total number of reads
plt.title(f'Base Frequency Logo (Vertically Stretched Alphabets)\nTotal Reads: {total_reads}', fontsize=14)
plt.ylim(0, 1.05)
plt.xlim(-0.5, num_plot_positions_range - 0.5)

base_font_size_pt = 100
color_map = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red', 'N': 'purple'}
base_y_position = 0.0

font_prop = FontProperties(family='sans-serif', weight='bold')

for x_idx, pos in enumerate(plot_positions_range):
    cumulative_y = base_y_position
    current_freqs = freqs_data.get(pos, {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0})
    sorted_bases = sorted(current_freqs.items(), key=lambda item: item[1], reverse=False)

    for base, freq in sorted_bases:
        if freq > 0:
            text_path = TextPath((0, 0), base, prop=font_prop, size=1)
            patch = PathPatch(text_path, facecolor=color_map[base], edgecolor='none', lw=0)

            bbox = text_path.get_extents()
            path_height = bbox.ymax - bbox.ymin
            scale_factor_x = letter_width / (bbox.xmax - bbox.xmin)
            scale_factor_y = freq / path_height if path_height != 0 else 1

            x_offset = x_idx - (letter_width / 2)

            transform = mpl.transforms.Affine2D().scale(scale_factor_x, scale_factor_y).translate(x_offset, cumulative_y)
            patch.set_transform(transform + ax.transData)

            ax.add_patch(patch)
            cumulative_y += freq

plt.tight_layout()

# Save the plot if an output file is specified
if output_file:
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
else:
    plt.show()
