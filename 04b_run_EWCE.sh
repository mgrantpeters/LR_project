#!/bin/bash

# Simple SLURM sbatch example
#SBATCH --job-name=EWCE
#SBATCH --ntasks=1
#SBATCH --time=3-00
#SBATCH --mem-per-cpu=40G
#SBATCH --partition=cpu

set -e -o pipefail

cd /nemo/lab/gandhis/home/users/grantpm/
ml purge 
ml Anaconda3
ml GEOS/3.9.1-GCC-11.2.0
. load_panpipes.sh

cd LR_project/
echo "Start running R"

Rscript 04b_EWCE_cell_specificity.R

echo "All done."