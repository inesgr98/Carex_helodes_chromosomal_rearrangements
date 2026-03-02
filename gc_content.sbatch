#!/bin/bash
#SBATCH --job-name=gc_content
#SBATCH --export=ALL
#SBATCH --output=gc_content_%j_output
#SBATCH --error=gc_content_%j_error
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00

### Set directories
INDIR=~/genome_assembly_helodes_thesis/C0_clean_assemblies


source ~/anaconda3/etc/profile.d/conda.sh
conda activate bedtools

#bedtools nuc -fi $INDIR/1JMC18.fasta -bed syn_100kb_window_1JMC18_affected_areas.bed > 100kb_window_1JMC18_affected_areas.gc.tsv
bedtools nuc -fi ~/genome_assembly_helodes_thesis/C1b_repeatmasker/4JMC19C_revcomp_mask.fasta -bed syn_100kb_window_4JMC19C_rearr_points_fixed.bed > 100kb_window_4JMC19C_rearr_points_fixed_mod.gc.tsv
conda deactivate

        

