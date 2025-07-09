import os
import glob
import numpy as np
import pandas as pd

# Set paths
INPUT_DIR = "/Users/921623492/Ecoli_Project/Data/merged_alignment"
OUTPUT_DIR = "/Users/921623492/Ecoli_Project/Data/MAF_output"
FASTA_PATH = "/Users/921623492/Ecoli_Project/Data/reduced_pangenome_blast.fa"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load all .npy paths
file_paths = glob.glob(os.path.join(INPUT_DIR, "*_pangenome_alignment.npy"))
if len(file_paths) == 0:
    raise FileNotFoundError("No .npy files found in the alignment directory!")

# Load gene names from FASTA
def extract_gene_names_from_fasta(fasta_path):
    with open(fasta_path, "r") as f:
        return [line[1:].strip() for line in f if line.startswith(">")]

gene_names = extract_gene_names_from_fasta(FASTA_PATH)[:15629]

# Combine .npy arrays into one matrix
combined_matrix = []
sample_ids = []

for path in file_paths:
    arr = np.load(path, allow_pickle=True)
    combined_matrix.append(arr)
    sample_id = os.path.basename(path).split(".")[0]
    sample_ids.append(sample_id)

# Create DataFrame
combined_df = pd.DataFrame(
    np.array(combined_matrix),
    index=sample_ids,
    columns=gene_names
)
combined_df.index.name = "Sample_ID"

# Save as CSV (or optionally Parquet for efficiency)
csv_path = os.path.join(OUTPUT_DIR, "combined_pangenome_matrix.csv")
combined_df.to_csv(csv_path)
print(f"Saved combined dataframe to: {csv_path}")
