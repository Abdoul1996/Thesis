#!/bin/bash
#SBATCH --job-name=combined_600_samples
#SBATCH --output=save_df_output.log
#SBATCH --error=save_df_error.log
#SBATCH --partition=cpucluster
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64


# Load environment if needed (uncomment and adjust)
# source ~/.bashrc
# conda activate abdoul_thesis

echo "Job started at $(date)"

python Script/save_combine.py

echo "Job completed at $(date)"