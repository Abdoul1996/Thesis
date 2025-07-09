import os
import glob
import numpy as np
import pandas as pd
from dask import delayed, compute
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster  # <-- NEW

# Start a local cluster using all 64 cores
cluster = LocalCluster(n_workers=64, threads_per_worker=1)
client = Client(cluster)
print(client)
print("Dask dashboard:", client.dashboard_link)


# Paths
INPUT_DIR = "/Users/921623492/Ecoli_Project/Data/merged_alignment"
OUTPUT_DIR = "/Users/921623492/Ecoli_Project/Experiment/results"
FASTA_PATH = "/Users/921623492/Ecoli_Project/Data/reduced_pangenome_blast.fa"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load gene names
def extract_gene_names_from_fasta(fasta_path):
    with open(fasta_path, "r") as f:
        return [line[1:].strip() for line in f if line.startswith(">")]

gene_names = extract_gene_names_from_fasta(FASTA_PATH)[:15629]

# Load all .npy paths
file_paths = glob.glob(os.path.join(INPUT_DIR, "*_pangenome_alignment.npy"))

# Delayed function to load a file and return a DataFrame
@delayed
def load_npy_to_df(path):
    arr = np.load(path, allow_pickle=True)
    sample_id = os.path.basename(path).split(".")[0]
    df = pd.DataFrame([arr], index=[sample_id], columns=gene_names)
    df.index.name = "Sample_ID"
    return df

# Wrap all file loads in Dask
delayed_dfs = [load_npy_to_df(path) for path in file_paths]

# Combine all with Dask
print("Launching parallel execution with Dask...")
combined_df = dd.from_delayed(delayed_dfs)

# Save to Parquet (fast + efficient)
output_path = os.path.join(OUTPUT_DIR, "combined_pangenome_matrix.parquet")
combined_df.to_parquet(output_path)
print(f"Dask output saved to: {output_path}")
