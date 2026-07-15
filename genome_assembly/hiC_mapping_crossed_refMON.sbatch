#!/bin/bash
#SBATCH --job-name=hiC_mapping_crossed_refMON
#SBATCH --output=hiC_mapping_crossed_refMON_%j_output
#SBATCH --error=hiC_mapping_crossed_refMON_%j_error
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
#SBATCH --time=23:30:00

#load directories
INDIR_hic=/home/igomez/genome_assembly_helodes_thesis/raw_data/hiC
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/hiC_mapping_crossed
SOFTWAREDIR=/home/igomez/software
REF=~/genome_assembly_helodes_thesis/C1b_repeatmasker/4JMC19C_revcomp_unmasked.fasta

##actívate conda
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ref_map
    ### Map HiC reads to ref using bwa mem

## First we need to map them as single reads, but we do it twice for both DNA strands

mkdir -p $OUTDIR/refMON/rawbamdir

bwa mem -t 16 $REF $INDIR_hic/1JMC18_R1_trimmed.fastq.gz |
samtools view -@ 16 -Sb - > $OUTDIR/refMON/rawbamdir/refMON_1.bam

bwa mem -t 16 $REF $INDIR_hic/1JMC18_R2_trimmed.fastq.gz |
samtools view -@ 16 -Sb - > $OUTDIR/refMON/rawbamdir/refMON_2.bam

echo "single read mapping done for refMON"

## We will eliminate the 3' end as it is redundant with the 5' end of the other strand.

mkdir -p $OUTDIR/refMON/filtbamdir

samtools view -h $OUTDIR/refMON/rawbamdir/refMON_1.bam | perl $SOFTWAREDIR/filter_five_end.pl | \
samtools view -Sb - > $OUTDIR/refMON/filtbamdir/refMON_1.bam
samtools view -h $OUTDIR/refMON/rawbamdir/refMON_2.bam | perl $SOFTWAREDIR/filter_five_end.pl | \
samtools view -Sb - > $OUTDIR/refMON/filtbamdir/refMON_2.bam

echo " 3' end eliminated from refMON "

## Now we will pair the two mapped DNA strands using  Map quality 0! is an important parameter which
##outputs a sorted, mapping quality filtered, paired-end BAM file. We then add read groups to this BAM
##file using Picard Tools.

    mkdir -p $OUTDIR/refMON/tmpdir

    perl $SOFTWAREDIR/two_read_bam_combiner.pl $OUTDIR/refMON/filtbamdir/refMON_1.bam \
    $OUTDIR/refMON/filtbamdir/refMON_2.bam samtools 1 | \
    samtools view -bS -t ${REF}.fai - | samtools sort -@ 16 -o $OUTDIR/refMON/tmpdir/refMON.bam -

    mkdir -p $OUTDIR/refMON/pairdir
    picard -Xmx100g AddOrReplaceReadGroups INPUT=$OUTDIR/refMON/tmpdir/refMON.bam OUTPUT=$OUTDIR/refMON/pairdir/refMON.bam ID=refMON LB=refMON SM=scaffolding PL=ILLUMINA PU=none

echo "read pairing done for refMON"

## Remove PCR duplicates from pair-end BAM file

    mkdir -p $OUTDIR/refMON/repdir
    picard -Xmx100g MarkDuplicates INPUT=$OUTDIR/refMON/pairdir/refMON.bam OUTPUT=$OUTDIR/refMON/repdir/refMON_final.bam METRICS_FILE=$OUTDIR/refMON/repdir/metrics.replicates.txt TMP_DIR=$OUTDIR/refMON/tmpdir ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

echo "PCR duplicates removed from refMON"

##Finally, index the final file and calculate statistics to monitor the mapping

samtools index $OUTDIR/refMON/repdir/refMON_final.bam

perl $SOFTWAREDIR/get_stats.pl $OUTDIR/refMON/repdir/refMON_final.bam > $OUTDIR/refMON/repdir/refMON_final.bam.stats

## Not only calculate statistics of the final file but also from intermediate steps. Add them to flagstat file
    touch $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats raw bam R1" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/rawbamdir/refMON_1.bam >> $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats raw bam R2" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/rawbamdir/refMON_2.bam >> $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats bam chimeric reads filtered R1" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/filtbamdir/refMON_1.bam >> $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats bam chimeric reads filtered R2" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/filtbamdir/refMON_2.bam >> $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats bam paired and mapping quality filtered" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/tmpdir/refMON.bam >> $OUTDIR/refMON/flagstat.txt
    echo ">>> Stats final bam (duplicates filtered)" >> $OUTDIR/refMON/flagstat.txt
    samtools flagstat $OUTDIR/refMON/repdir/refMON_final.bam >> $OUTDIR/refMON/flagstat.txt

echo "index and stats done for refMON"


conda deactivate

### SOFTWARE VERSIONS
# fastp 0.23.4
# bwa 0.7.15
# Picard 1.141
# Arima mapping pipeline V03
