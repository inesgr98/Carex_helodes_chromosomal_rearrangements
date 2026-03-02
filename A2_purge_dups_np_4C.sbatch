#!/usr/bin/bash
#SBATCH --job-name=A2_purge_dups_after_hifiasm_4C
#SBATCH --output=A2_purge_dups_after_hifiasm_4C_output_%j
#SBATCH --error=A2_purge_dups_after_hifiasm_4C_error_%j
#SBATCH --nodes=1                 
#SBATCH --ntasks=48
#SBATCH --mem=60G
#SBATCH --partition=standard
#SBATCH --time=23:30:00

# Activate the conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate purge_duplicates

# Define variables
### Create function with all the purge_dups commands
	REF=/home/igomez/genome_assembly_helodes_thesis/A1_hifiasm/4JMC19C/4JMC19C.hifihic.asm.hic.p_ctg.fa
	ALT=/home/igomez/genome_assembly_helodes_thesis/A1_hifiasm/4JMC19C/4JMC19C.hifihic.asm.hic.a_ctg.fa
	READS=/home/igomez/genome_assembly_helodes_thesis/raw_data/PacBio/4JMC19C.fastq
	OUTDIR="/home/igomez/genome_assembly_helodes_thesis/A2_purge_dups"

# Step1: Map reads to reference and calculate coverage
	mkdir -p $OUTDIR/4JMC19C
	mkdir -p $OUTDIR/4JMC19C
	cd $OUTDIR/4JMC19C
	#minimap2 -xasm20 -t 12 $REF $READS | gzip -c - > $OUTDIR/4JMC19C/4JMC19C.paf.gz
	#pbcstat $OUTDIR/4JMC19C/4JMC19C.paf.gz
	#calcuts $OUTDIR/4JMC19C/PB.stat > $OUTDIR/4JMC19C/cutoffs \
	#2>$OUTDIR/4JMC19C/calcults.log
	#echo "step 1 4JMC19C done"

# Generate histogram
    	
	#hist_plot.py -c cutoffs PB.stat cov_unpurged_hist.png
	
# Step2: Split the assembly to do self-self alignment
	#split_fa $REF > $OUTDIR/4JMC19C/4JMC19C.split
	#minimap2 -xasm5 -t 12 -DP $OUTDIR/4JMC19C/4JMC19C.split $OUTDIR/4JMC19C/4JMC19C.split | gzip -c - \
	#> $OUTDIR/4JMC19C/4JMC19C.split.self.paf.gz
	#echo "step 2 4JMC19C done"

# Step3: purge haplotigs and overlaps
	#purge_dups -2 -T $OUTDIR/4JMC19C/cutoffs -c $OUTDIR/4JMC19C/PB.base.cov $OUTDIR/4JMC19C/4JMC19C.split.self.paf.gz \
	#> $OUTDIR/4JMC19C/dups.bed 2> $OUTDIR/4JMC19C/purge_dups.log
	echo "step 3 4JMC19C done"

# Step4: get purged primary and haplotigs sequences
	#get_seqs -e $OUTDIR/4JMC19C/dups.bed $REF
	echo "step 4 4JMC19C done"

# Step5: Merge alternate haplotype and hap.fa to get a decent haplotig set and repeat steps above on ALTREF
	mkdir -p $OUTDIR/4JMC19C/alternate
	cd $OUTDIR/4JMC19C/alternate
	#cat $ALT $OUTDIR/4JMC19C/hap.fa > $OUTDIR/4JMC19C/alternate/4JMC19C_alt_hap.fa
	ALTREF=$OUTDIR/4JMC19C/alternate/4JMC19C_alt_hap.fa
	minimap2 -xasm20 -t 12 $ALTREF $READS | gzip -c - > $OUTDIR/4JMC19C/alternate/4JMC19C.paf.gz
	pbcstat $OUTDIR/4JMC19C/alternate/4JMC19C.paf.gz
	calcuts $OUTDIR/4JMC19C/alternate/PB.stat > $OUTDIR/4JMC19C/alternate/cutoffs \
	2>$OUTDIR/4JMC19C/alternate/calcults.log
echo "step 1 part 2 4JMC19C done"
	split_fa $ALTREF > $OUTDIR/4JMC19C/alternate/4JMC19C.split
	minimap2 -xasm5 -t 12 -DP $OUTDIR/4JMC19C/alternate/4JMC19C.split \
	$OUTDIR/4JMC19C/alternate/4JMC19C.split | gzip -c - \
	> $OUTDIR/4JMC19C/alternate/4JMC19C.split.self.paf.gz
echo "step 2 part 2 4JMC19C done"
purge_dups -2 -T $OUTDIR/4JMC19C/alternate/cutoffs -c $OUTDIR/4JMC19C/alternate/PB.base.cov $OUTDIR/4JMC19C/alternate/4JMC19C.split.self.paf.gz > $OUTDIR/4JMC19C/alternate/dups.bed 2> $OUTDIR/4JMC19C/alternate/purge_dups.log
echo "step 3 part 2 4JMC19C done"

	get_seqs -e $OUTDIR/4JMC19C/alternate/dups.bed $ALTREF
echo "step 4 part 2 4JMC19C done"
	


### SOFTWARE VERSIONS
# minimap2 2.28
# purge_dups 1.2.6



