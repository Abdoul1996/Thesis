import os
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm

# ---------------------
# Helper Functions
# ---------------------
def extract_gene_names_from_fasta(fasta_path):
    with open(fasta_path, "r") as f:
        return [line[1:].strip() for line in f if line.startswith(">")]

# ---------------------
# Main Logic
# ---------------------
if __name__ == "__main__":

    # Paths
    INPUT_DIR = "/Users/921623492/Ecoli_Project/Data/merged_alignment"
    OUTPUT_DIR = "/Users/921623492/Ecoli_Project/Experiment/results"
    FASTA_PATH = "/Users/921623492/Ecoli_Project/Data/reduced_pangenome_blast.fa"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("📂 Loading gene names and file paths...")
    gene_names = extract_gene_names_from_fasta(FASTA_PATH)[:15629]
    file_paths = glob.glob(os.path.join(INPUT_DIR, "*_pangenome_alignment.npy"))
    print(f"🧬 Found {len(file_paths)} .npy files")

    combined_matrix = []
    sample_ids = []

    # Process files with progress bar
    for path in tqdm(file_paths, desc="🔄 Processing alignment files"):
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

    # Save
    output_path = os.path.join(OUTPUT_DIR, "combined_pangenome_matrix.csv")
    combined_df.to_csv(output_path)
    print(f"Done! File saved at: {output_path}")
