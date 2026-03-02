#!/bin/bash
#SBATCH --job-name=B1_hiC_trimming_np
#SBATCH --output=B1_hiC_trimming_output_np_%j
#SBATCH --error=B1_hiC_trimming_error_np_%j
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
#SBATCH --time=23:30:00

####Set directories
INDIR_hic=/home/igomez/genome_assembly_helodes_thesis/raw_data/hiC
SOFTWAREDIR=/home/igomez/software
cd $INDIR_hic

##Set conda source and activate environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ref_map

## Trimming 5 bases out of 5' end out of read 1 and read 2,
##as they contain adapter sequencuences
gunzip -c $INDIR_hic/1JMC18_hiC_1.fastq.gz | \
awk '{ if(NR%2==0) {print substr($1,6)} else {print} }' | \
gzip > $INDIR_hic/1JMC18_R1_5trim.fastq.gz

gunzip -c $INDIR_hic/1JMC18_hiC_2.fastq.gz | \
awk '{ if(NR%2==0) {print substr($1,6)} else {print} }' | \
gzip > $INDIR_hic/1JMC18_R2_5trim.fastq.gz 

echo "5 bases trim done for 1JMC18"

##Trim poly-G tails from the ends of sequencing reads,as they are sequencing artifacts
    
fastp --trim_poly_g -l 120 -q 30 -w 16 \
--in1 $INDIR_hic/1JMC18_R1_5trim.fastq.gz \
--out1 $INDIR_hic/1JMC18_R1_trimmed.fastq.gz \
--in2 $INDIR_hic/1JMC18_R2_5trim.fastq.gz \
--out2 $INDIR_hic/1JMC18_R2_trimmed.fastq.gz

echo "poly-G tail trim done for 1JMC18"

conda deactivate

