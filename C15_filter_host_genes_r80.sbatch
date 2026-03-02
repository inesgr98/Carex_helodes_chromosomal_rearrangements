#!/bin/bash
#SBATCH --job-name=C15_blast_host_genes
#SBATCH --output=C15_blast_host_genes_%j.out
#SBATCH --error=C15_blast_host_genes_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=80G
#SBATCH --time=48:00:00


### Check if sequences from the EarlGrey library are in fact host genes
### By blasting the library against the nr Arthropoda database
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes
LIB=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes/C_helodes_lib_r_clstrd80.fasta ### LIB from EarlGrey output

###Actiavte conda
source ~/anaconda3/etc/profile.d/conda.sh
#conda activate diamond 

cd $OUTDIR

#diamond makedb -p 48 --in Poales_nr.fasta -d Poales_nr_diamond.db
#diamond blastx -p 48 -d Poales_nr_diamond.db -q $LIB -o C_helodes_lib_blast.tsv

### Filter the blast hits so that we keep only one hit per repeat
#awk -F'\t' '!seen[$1]++' C_helodes_lib_blast.tsv > C_helodes_lib_blast_uniq.tsv

### Get the list of GenBank IDs corresponding to the hits
#cat C_helodes_lib_blast_uniq.tsv | cut -f 2 > C_helodes_lib_blast_uniq_GenBankID.tsv
#conda deactivate

#echo "GenBank ID list generated"


### Conver the list of GenBank IDs to protein fasta sequences

conda activate entrez-direct
#while IFS= read -r line <&3; do
#    esearch -db protein -query "$line" | efetch -format fasta >> C_helodes_lib_blast_uniq_GenBank.fasta
#done 3< C_helodes_lib_blast_uniq_GenBankID.tsv
### Too long for batch CD search, separate into 6
awk -v n=$(grep -c '^>' C_helodes_lib_blast_uniq_GenBank.fasta) -v f=C_helodes_lib_blast_uniq_GenBank.fasta \
    'BEGIN { p = int(n / 8); r = n % 8; split_count = p + (r > 0 ? 1 : 0); extra = r - (r > 0 ? 1 : 0); part = 1; seq = 0; file = "lib_part_" part ".fasta" } /^>/ { if (seq >= split_count) { part++; seq = 0; split_count = p + (part <= r ? 1 : 0); file = "lib_part_" part ".fasta" } seq++ } { print >> file }' C_helodes_lib_blast_uniq_GenBank.fasta
conda deactivate


### Then CD search online on https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
### Download the results and look if any match to the known TE domains from Wicker et al.
### Delete the lines where there is a known TE domain to have the final list of GenBank IDs corresponding to host genes
### The list of blast matches that are probably host genes is lib_filtered_host_genes.tsv
### Now keep only the blast matches in the file C_helodes_lib_blast_uniq.tsv that really correspond to host genes


