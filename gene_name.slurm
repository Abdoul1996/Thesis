#!/bin/bash
#SBATCH --job-name=blast-one          # Optional: name your job
#SBATCH --nodes=1                     # Number of nodes (no space after '=')
#SBATCH --partition=cputest             # Your target partition (e.g., cpu, general)
#SBATCH --output=/Users/921623492/Ecoli_Project/Experiment/logs/blast_one_%j.out
#SBATCH --error=/Users/921623492/Ecoli_Project/Experiment/logs/blast_one_%j.err
# Load Conda (if your HPC requires it — skip if already available)

module load anaconda || source ~/miniconda3/etc/profile.d/conda.sh

# Activate your environment
conda activate abdoul_thesis

python gene_names.py                # Run your script