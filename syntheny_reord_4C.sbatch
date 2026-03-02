#!/bin/bash
#SBATCH --job-name=syntheny_reord_scaffold_4C
#SBATCH --output=syntheny_reord_scaffold_4C_%j.out
#SBATCH --error=syntheny_reord_scaffold_4C_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=6:00:00

source ~/anaconda3/etc/profile.d/conda.sh
conda activate seqtik_env

input_fasta=~/genome_assembly_helodes_thesis/C1b_repeatmasker/4JMC19C/4JMC19C.fasta.masked
output_fasta=~/genome_assembly_helodes_thesis/C1b_repeatmasker/4JMC19C_revcomp_mask.fasta

##Name chromosomes to revert
#chromosomes="scaffold_24 scaffold_4 scaffold_7 scaffold_8 scaffold_13 scaffold_10 scaffold_14 scaffold_20 scaffold_23 scaffold_32 scaffold_30 scaffold_26 scaffold_27 scaffold_35 scaffold_33 scaffold_16 scaffold_36"

# Step 1: Generate chromosome list
#printf "scaffold_24\nscaffold_4\nscaffold_7\nscaffold_8\nscaffold_13\nscaffold_10\nscaffold_14\nscaffold_20\nscaffold_23\nscaffold_32\nscaffold_30\nscaffold_26\nscaffold_27\nscaffold_35\nscaffold_33\nscaffold_16\nscaffold_36" > chrom_list.txt

# Step 2: Reverse complement only the selected chromosomes
#seqkit grep -n -f chrom_list.txt --id-regexp "^(scaffold_[0-9]+)$" $input_fasta | seqkit seq -r -p -t DNA > reversed_chroms.fasta

# Step 3: Extract unmodified chromosomes
#seqkit grep -v -n -f chrom_list.txt $input_fasta > unchanged_chroms.fasta

# Step 4: Merge reversed and unchanged sequences
#cat reversed_chroms.fasta unchanged_chroms.fasta > $output_fasta

#comprobar archivos
seqkit stat $output_fasta
seqkit stat $input_fasta
conda deactivate

echo "Reverse complementation complete. Output saved in modified_genome.fasta"

