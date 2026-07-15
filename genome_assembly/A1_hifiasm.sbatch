#!/usr/bin/bash
#SBATCH --job-name=A1_hifiasm_hifihic
#SBATCH --output=A1_hifiasm_hifihic_output_%j
#SBATCH --error=A1_hifiasm_hifihic_error_%j
#SBATCH --nodes=1                  
#SBATCH --ntasks=32
#SBATCH --mem=40G
#SBATCH --partition=standard
#SBATCH --time=12:00:00


# Define variables
INDIR_hifi="/home/igomez/genome_assembly_helodes_thesis/raw_data/PacBio"
INDIR_hic="/home/igomez/genome_assembly_helodes_thesis/raw_data/hiC"
OUTDIR="/home/igomez/genome_assembly_helodes_thesis/A1_hifiasm"
SAMPLE_DIR="/home/igomez/genome_assembly_helodes_thesis/parents.txt"

#Run hifiasm
hifiasm_assembly() {
SAMPLE=$1
mkdir -p $OUTDIR/${SAMPLE}
# Activate the conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate hifiasm
hifiasm -t 32 -D 10 -N 200 -o $OUTDIR/${SAMPLE}/${SAMPLE}.hifihic.asm --h1 $INDIR_hic/${SAMPLE}_hiC_1.fastq.gz --h2 $INDIR_hic/${SAMPLE}_hiC_2.fastq.gz $INDIR_hifi/${SAMPLE}.fastq
echo "hifiasm ${SAMPLE} done"
conda deactivate
}


### Parallelize over samples
export INDIR_hifi INDIR_hic OUTDIR SAMPLE_DIR
export -f hifiasm_assembly
parallel 'hifiasm_assembly {1}' :::: $SAMPLE_DIR

###   VERSIONS
# hifiasm 0.19.9
