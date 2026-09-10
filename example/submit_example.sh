#!/bin/bash -l
#SBATCH -p standard   # Partition name 
#SBATCH -n 1     # Number of cores
#SBATCH -J mixkin # Job name
#SBATCH --output=slurm-%j-%N-%u.txt
#SBATCH --error=slurm-%j-%N-%u.err

source mixkin_env/bin/activate
# conda activate mixkin_env

mixkin \
  -i example \
  -q example.Q \
  -y 2 \
  -o kinship_results.tsv

