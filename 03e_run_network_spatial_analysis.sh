#!/bin/bash

# Simple SLURM sbatch example
#SBATCH --job-name=LR_networks
#SBATCH --ntasks=1
#SBATCH --time=3-00
#SBATCH --mem-per-cpu=48G
#SBATCH --partition=cpu

set -e -o pipefail

cd /nemo/lab/gandhis/home/users/grantpm/
ml purge 
. load_panpipes.sh

cd LR_project/
echo "Start running python"

python3 03e_network_spatial_analysis.py

echo "All done."