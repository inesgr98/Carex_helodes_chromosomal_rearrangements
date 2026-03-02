#!/bin/bash
#SBATCH --job-name=B11_map_hiC_reads_step1
#SBATCH --output=B11_map_hiC_reads_step1_%j_output
#SBATCH --error=B11_map_hiC_reads_step1_%j_error
#SBATCH --export=ALL
#SBATCH --nodes=1                
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard
#SBATCH --time=12:00:00

##actívate conda 
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ref_map
SAMPLE_DIR="/home/igomez/genome_assembly_helodes_thesis/parents.txt"

mapHiCreads_step1() {
SAMPLE=$2

    ### Index ref, both with samtools and bwa
    REF=/home/igomez/genome_assembly_helodes_thesis/A2_purge_dups/${SAMPLE}/purged.fa
    bwa index -a bwtsw $REF
    samtools faidx $REF
}

### Parallelize over individuals
export SAMPLE_DIR
export -f mapHiCreads_step1
parallel 'mapHiCreads_step1 {1}' :::: $SAMPLE_DIR

### SOFTWARE VERSIONS
# fastp 0.23.4
# bwa 0.7.15
# Picard 1.141
# Arima mapping pipeline V03




