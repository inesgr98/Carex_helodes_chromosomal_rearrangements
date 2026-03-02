#!/usr/bin/bash
#SBATCH --job-name=A2_purge_dups_after_hifiasm
#SBATCH --output=A2_purge_dups_after_hifiasm_output_%j
#SBATCH --error=A2_purge_dups_after_hifiasm_error_%j
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
	REF=/home/igomez/genome_assembly_helodes_thesis/A1_hifiasm/1JMC18/1JMC18.hifihic.asm.hic.p_ctg.fa
	ALT=/home/igomez/genome_assembly_helodes_thesis/A1_hifiasm/1JMC18/1JMC18.hifihic.asm.hic.a_ctg.fa
	READS=/home/igomez/genome_assembly_helodes_thesis/raw_data/PacBio/1JMC18.fastq
	OUTDIR="/home/igomez/genome_assembly_helodes_thesis/A2_purge_dups"

# Step1: Map reads to reference and calculate coverage
	mkdir -p $OUTDIR/1JMC18
	mkdir -p $OUTDIR/1JMC18
	cd $OUTDIR/1JMC18
	minimap2 -xasm20 -t 12 $REF $READS | gzip -c - > $OUTDIR/1JMC18/1JMC18.paf.gz
	pbcstat $OUTDIR/1JMC18/1JMC18.paf.gz
	calcuts $OUTDIR/1JMC18/PB.stat > $OUTDIR/1JMC18/cutoffs \
	2>$OUTDIR/1JMC18/calcults.log
	echo "step 1 1JMC18 done"

# Step2: Split the assembly to do self-self alignment
	split_fa $REF > $OUTDIR/1JMC18/1JMC18.split
	minimap2 -xasm5 -t 12 -DP $OUTDIR/1JMC18/1JMC18.split $OUTDIR/1JMC18/1JMC18.split | gzip -c - \
	> $OUTDIR/1JMC18/1JMC18.split.self.paf.gz
	echo "step 2 1JMC18 done"

# Step3: purge haplotigs and overlaps
	purge_dups -2 -T $OUTDIR/1JMC18/cutoffs -c $OUTDIR/1JMC18/PB.base.cov $OUTDIR/1JMC18/1JMC18.split.self.paf.gz \
	> $OUTDIR/1JMC18/dups.bed 2> $OUTDIR/1JMC18/purge_dups.log
	echo "step 3 1JMC18 done"

# Step4: get purged primary and haplotigs sequences
	get_seqs -e $OUTDIR/1JMC18/dups.bed $REF
	echo "step 4 1JMC18 done"

# Step5: Merge alternate haplotype and hap.fa to get a decent haplotig set and repeat steps above on ALTREF
	mkdir -p $OUTDIR/1JMC18/alternate
	cd $OUTDIR/1JMC18/alternate
	cat $ALT $OUTDIR/1JMC18/hap.fa > $OUTDIR/1JMC18/alternate/1JMC18_alt_hap.fa
	ALTREF=$OUTDIR/1JMC18/alternate/1JMC18_alt_hap.fa
	minimap2 -xasm20 -t 12 $ALTREF $READS | gzip -c - > $OUTDIR/1JMC18/alternate/1JMC18.paf.gz
	pbcstat $OUTDIR/1JMC18/alternate/1JMC18.paf.gz
	calcuts $OUTDIR/1JMC18/alternate/PB.stat > $OUTDIR/1JMC18/alternate/cutoffs \
	2>$OUTDIR/1JMC18/alternate/calcults.log
echo "step 1 part 2 1JMC18 done"
	split_fa $ALTREF > $OUTDIR/1JMC18/alternate/1JMC18.split
	minimap2 -xasm5 -t 12 -DP $OUTDIR/1JMC18/alternate/1JMC18.split \
	$OUTDIR/1JMC18/alternate/1JMC18.split | gzip -c - \
	> $OUTDIR/1JMC18/alternate/1JMC18.split.self.paf.gz
echo "step 2 part 2 1JMC18 done"
purge_dups -2 -T $OUTDIR/1JMC18/alternate/cutoffs -c $OUTDIR/1JMC18/alternate/PB.base.cov $OUTDIR/1JMC18/alternate/1JMC18.split.self.paf.gz > $OUTDIR/1JMC18/alternate/dups.bed 2> $OUTDIR/1JMC18/alternate/purge_dups.log
echo "step 3 part 2 1JMC18 done"

	get_seqs -e $OUTDIR/1JMC18/alternate/dups.bed $ALTREF
echo "step 4 part 2 1JMC18 done"
	


### SOFTWARE VERSIONS
# minimap2 2.28
# purge_dups 1.2.6



