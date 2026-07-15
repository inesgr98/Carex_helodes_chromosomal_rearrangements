#!/usr/bin/bash
#SBATCH --job-name=unmask_4JMC19C
#SBATCH --output=unmask_4JMC19C_%j_output
#SBATCH --error=unmask_4JMC19C_%j_error
#SBATCH --nodes=1                  # Using one node
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --partition=standard
#SBATCH --time=12:00:00

source ~/anaconda3/etc/profile.d/conda.sh
#conda activate seqtik_env
#seqkit seq -u 4JMC19C_revcomp_mask.fasta > 4JMC19C_revcomp_unmasked.fasta
#conda deactivate

conda activate ref_map
#samtools faidx 4JMC19C_revcomp_unmasked.fasta

bwa index 4JMC19C_revcomp_unmasked.fasta
conda deactivate
