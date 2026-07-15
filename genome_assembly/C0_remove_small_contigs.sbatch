#!/bin/bash
#SBATCH --job-name=C0_clean_assemblies
#SBATCH --export=ALL
#SBATCH --output=C0_clean_assemblies_%j_output
#SBATCH --error=C0_clean_assemblies_%j_error
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00

### Set directories
INDIR=~/genome_assembly_helodes_thesis/B4_decontaminate
OUTDIR=~/genome_assembly_helodes_thesis/C0_clean_assemblies
SAMPLE_DIR=~/genome_assembly_helodes_thesis/parents.txt



for SAMPLE in $(cat $SAMPLE_DIR | cut -f1)
do


source ~/anaconda3/etc/profile.d/conda.sh
conda activate seqtik_env

MINLEN=4000000 #maybe modify it according to contig length list, hiC and telomeres. 
seqkit seq --threads 16 --min-len $MINLEN $INDIR/${SAMPLE}/eukarya_${SAMPLE}_JBAT.FINAL.fa > $OUTDIR/${SAMPLE}.fasta
conda deactivate

echo "clean assemblies generated without small, umapped scaffolds for ${SAMPLE}"
        
done

