#!/bin/bash
#SBATCH --job-name=syntheny_helodes
#SBATCH --output=syntheny_helodes_%j.out
#SBATCH --error=syntheny_helodes_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --time=48:00:00


### Check if sequences from the EarlGrey library are in fact host genes
### By blasting the library against the nr Arthropoda database
INDIR=/home/igomez/genome_assembly_helodes_thesis/C2_braker
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/syntheny

###Actiavte conda
source ~/anaconda3/etc/profile.d/conda.sh

conda activate r-env
Rscript test_packages.R
conda deactivate


