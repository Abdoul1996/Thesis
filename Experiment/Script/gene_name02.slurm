#!/bin/bash
#SBATCH --job-name=blast_batch
#SBATCH --partition=cputest 
#SBATCH --output=/Users/921623492/Ecoli_Project/Experiment/logs/blast_batch_%j.out
#SBATCH --error=/Users/921623492/Ecoli_Project/Experiment/logs/blast_batch_%j.err
#SBATCH --ntasks=1

# Run the Python script
python /Users/921623492/Ecoli_Project/gene_name02.py \
  --input_csv /Users/921623492/Ecoli_Project/Experiment/results/combined_pangenome_matrix.csv \
  --genes FAHFDEJI_02221 BPDJECED_04485 GMELEDML_04617 \
  --output_csv /Users/921623492/Ecoli_Project/Experiment/results/batch_results_1.csv
