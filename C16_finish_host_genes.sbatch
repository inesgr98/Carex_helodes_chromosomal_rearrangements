#!/bin/bash
#SBATCH --job-name=C16_finish_host_genes
#SBATCH --output=C16_finish_host_genes_%j.out
#SBATCH --error=C16_finish_host_genes_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00


### Check if sequences from the EarlGrey library are in fact host genes
### By blasting the library against the nr Arthropoda database
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes
LIB=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes/C_helodes_lib_r_clstrd80.fasta ### LIB from EarlGrey output

###Actiavte conda
source ~/anaconda3/etc/profile.d/conda.sh
 

cd $OUTDIR

### Then CD search online on https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
### Download the results and look if any match to the known TE domains from Wicker et al.
### Delete the lines where there is a known TE domain to have the final list of GenBank IDs corresponding to host genes
### The list of blast matches that are probably host genes is lib_filtered_host_genes.tsv
### Now keep only the blast matches in the file C_helodes_lib_blast_uniq.tsv that really correspond to host genes
awk 'FNR==NR { keep[$1]; next } ($2 in keep)' lib_corr_id.tsv C_helodes_lib_blast_uniq.tsv > C_helodes_lib_blast_only_hostgenes.tsv
### Now remove the host genes from the Erebia lib (keep in $LIB only the sequences that have a header that is not present in the first column of Erebia_lib_blast_filtered_only_host_gene.tsv)
awk 'BEGIN { FS="\t" } 
    FNR==NR { exclude[$1]; next } 
    /^>/ { 
        header=$0; 
        sub(/^>/, "", header);  # Remove leading ">" from the header
        gsub(/[ \t]+$/, "", header);  # Trim trailing whitespace
        keep = !(header in exclude);
    } 
    keep { print $0 }'  lib_part_corr.fasta $LIB > C_helodes_lib_clstrd80_hostgenesfilt.fasta


