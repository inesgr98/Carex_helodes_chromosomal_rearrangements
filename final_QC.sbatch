#!/bin/bash
#SBATCH --job-name=final_QC_%j
#SBATCH --export=ALL
#SBATCH --output=final_QC_%j_output
#SBATCH --error=final_QC_%j_error
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00

#set source for conda
source ~/anaconda3/etc/profile.d/conda.sh

### Set directories
OUTDIR="/home/igomez/genome_assembly_helodes_thesis/C0_clean_assemblies/final_QC"
INDIR="/home/igomez/genome_assembly_helodes_thesis/C0_clean_assemblies"
SAMPLE_DIR="/home/igomez/genome_assembly_helodes_thesis/parents.txt"

for SAMPLE in $(cat $SAMPLE_DIR | cut -f1)
do

REF=$INDIR/${SAMPLE}.fasta
mkdir -p $OUTDIR/${SAMPLE}
cd $OUTDIR/${SAMPLE}

## Find telomeres
conda activate tidk
tidk find -c Poales -d $OUTDIR/${SAMPLE} -o ${SAMPLE}_telomeres --log $REF
tidk plot -o ${SAMPLE}_telo_plot --tsv ${SAMPLE}.tsv
conda deactivate
echo "telomeric analysis done for ${SAMPLE}"

## Test completedness with Compleasm 
conda activate compleasm
compleasm run -a $REF -o $OUTDIR/${SAMPLE} -l poales_odb10 -t 16
conda deactivate
echo "compleasm done for ${SAMPLE}"

## Test contiguity with Gfastats
conda activate gfastats
gfastats $REF -j 16 > $OUTDIR/${SAMPLE}/${SAMPLE}_gfastats.txt
conda deactivate
echo "gfastats done for ${SAMPLE}"

## Get a list of contigs and scaffolds with their length
perl /home/igomez/software/contig_length_fasta.pl $REF > $OUTDIR/${SAMPLE}/${SAMPLE}_contig_length.txt
echo "contig list complete for ${SAMPLE}"

## Repetedness with merqury
#create directory to run merqury
mkdir -p $OUTDIR/${SAMPLE}/merqury
cd $OUTDIR/${SAMPLE}/merqury

#The meryl reads is calculated from the PacBio reads so it is the same for all the merqury analysis
#in the assembly
cp -r ~/genome_assembly_helodes_thesis/A3_QC/${SAMPLE}/${SAMPLE}.meryl .

# Set path to to the alternate assembly
REF_ALT=~/genome_assembly_helodes_thesis/A2_purge_dups/${SAMPLE}/alternate/purged_alt.fa

#Run merqury
conda activate merqury
merqury.sh ../${SAMPLE}.meryl $REF $REF_ALT final_merqury_${SAMPLE}
conda deactivate

done

### SOFTWARE VERSIONS
# compleasm 0.2.6
# contig_length_fasta.pl https://github.com/ishengtsai/perlscripts_jit/blob/master/contig_length_fasta.pl
# gfastat 1.3.6
# merqury 1.3
# tidk 0.2.41

