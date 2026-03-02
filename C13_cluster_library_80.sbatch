#!/bin/bash
#SBATCH --job-name=C13_cluster_repeats
#SBATCH --output=C13_cluster_repeats_%j.out
#SBATCH --error=C13_cluster_repeats_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=100:00:00

OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes
INDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library ### LIB from EarlGrey output

###Actiavte conda
source ~/anaconda3/etc/profile.d/conda.sh
conda activate cd_hit

cd-hit-est -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -r 1 -T 24 -M 0 \
-i $INDIR/1JMC18/1JMC18r1_conc_EarlGrey/1JMC18r1_conc_Curated_Library/1JMC18r1_conc_combined_library.fasta -o $OUTDIR/C_helodes_lib_r_clstrd80.fasta
### Using a 95% sequence identity cut-of f leads to reducing the library from 14308 to 11558 sequences
### Still way too much, use the standard 80% at the risk of blurring some sub families
### Using 80% reduces the library to 9864 sequences
conda deactivate

