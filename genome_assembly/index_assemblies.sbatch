#!/bin/bash
#SBATCH --job-name=C0_index_assemblies
#SBATCH --export=ALL
#SBATCH --output=C0_index_assemblies_%j_output
#SBATCH --error=C0_index_assemblies_%j_error
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

### Set directorie

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ref_map

  # REF=/home/igomez/holocentromere_repeats/index_pilulifera/GCA_964058765.1_lpCarAren1.1_genomic.fna

#    bwa index -a bwtsw $REF
    samtools faidx 1JMC18/1JMC18.fasta.masked
    conda deactivate    


