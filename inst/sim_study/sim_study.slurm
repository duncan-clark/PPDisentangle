#!/bin/bash
#SBATCH --job-name=PPDisentangle_sim
#SBATCH --output=sim_study_%j.out
#SBATCH --error=sim_study_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

module load R 2>/dev/null || true

Rscript sim_study.R
