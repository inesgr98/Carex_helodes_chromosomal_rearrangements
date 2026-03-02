#!/usr/bin/bash
#SBATCH --job-name=bed_extract_4JMC19C
#SBATCH --output=bed_extract_reord_%j_output
#SBATCH --error=bed_extract_reord_%j_error
#SBATCH --nodes=1                  # Using one node
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --partition=standard
#SBATCH --time=12:00:00

source ~/anaconda3/etc/profile.d/conda.sh
conda activate bedtools

FASTA=4JMC19C_revcomp_mask.fasta
OUTBED=repeats_4JMC19C_revcomp.bed

awk '
BEGIN{OFS="\t"}
/^>/{
  chr=substr($0,2)
  pos=0
  next
}
{
  line=$0
  while (match(line, /[acgt]+/)) {
    start = pos + RSTART - 1
    end   = start + RLENGTH
    print chr, start, end
    line = substr(line, RSTART + RLENGTH)
    pos += RSTART + RLENGTH - 1
  }
  pos += length($0)
}
' $FASTA > $OUTBED

sort -k1,1 -k2,2n $OUTBED > repeats_4JMC19C_revcomp.sorted.bed

bedtools merge -i repeats_4JMC19C_revcomp.sorted.bed > repeats_4JMC19C_revcomp.merged.bed

conda deactivate
