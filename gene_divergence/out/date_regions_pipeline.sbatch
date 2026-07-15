#!/bin/bash
#SBATCH --job-name=divergence_translocation
#SBATCH --output=divergence_translocation_%j.log
#SBATCH --error=divergence_translocation_%j.err
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH --mem=16G

source ~/anaconda3/etc/profile.d/conda.sh

###Extract CDS and genes
#input#
GFF3_FILE="/home/igomez/genome_assembly_helodes_thesis/C2b_braker/4JMC19C/braker/braker.gff3"
SYNTENY_HITS="/home/igomez/genome_assembly_helodes_thesis/C3b_syntheny/syntenicHits/C_helodes_4JMC19C_vs_C_helodes_1JMC18.synHits.txt"
GENOME_4JMC19C_FA="/home/igomez/genome_assembly_helodes_thesis/C1b_repeatmasker/4JMC19C_revcomp_mask.fasta"
GENOME_1JMC18_FA="/home/igomez/genome_assembly_helodes_thesis/C1b_repeatmasker/1JMC18/1JMC18.fasta.masked"
GFF3_FILE_1JMC18="/home/igomez/genome_assembly_helodes_thesis/C2b_braker/1JMC18/braker/braker.gff3"
GENE_LIST="target_genes_4C.txt"
CDS_LIST="target_cds_4C.txt"
CHR_NAME="scaffold_24"
MAP_DELIMITER='\t'

# Inversion coordinates for Carex borbonica

START_COORD="2986606"
END_COORD="3679665"

echo "--- Starting Gene ID Extraction (Outside Inversion) ---"
echo "Targeting features on Scaffold: ${CHR_NAME} *INSIDE* the range [${START_COORD} - ${END_COORD}]."

# Use awk with a forced tab delimiter (-F'\t') for GFF3 standard compliance.
awk -F'\t' -v chr="$CHR_NAME" -v start="$START_COORD" -v end="$END_COORD" '
# 1. Filter by Chromosome Name
$1 == chr &&
# 2. Filter by feature type "gene"
$3 == "gene" &&
# 3. CRITICAL: Filter by coordinates (start and end must be INSIDE the inversion boundaries)
# A gene is inside if its END ($5) is before the inversion START ($4 < start)
# OR if its START ($4) is after the inversion END ($4 > end).
( $5 < start || $4 > end ) {

    # 4. Extract and clean the Gene ID from the 9th column (Attributes)
    if (match($9, /ID=([^;]+)/)) {
        # RSTART is the start position of the match, RLENGTH is the length.
        # We start 3 characters after RSTART to skip the "ID=" prefix
        gene_id = substr($9, RSTART + 3, RLENGTH - 3)
        print gene_id
    }
}' "$GFF3_FILE" > "$GENE_LIST"

# --- Report Results ---
COUNT=$(wc -l < "$GENE_LIST")

if [ $? -eq 0 ]; then
    echo "Successfully extracted ${COUNT} gene IDs inside the inversion."
    echo "Gene list saved to: ${GENE_LIST}"
    echo "--- Extraction Complete ---"
else
    echo "ERROR: Awk command failed. Check your GFF3 file format."
fi


echo "--- Starting CDS ID Extraction (Outside Inversion) ---"
echo "Targeting features on Scaffold: ${CHR_NAME} *INSIDE* the range [${START_COORD} - ${END_COORD}]."

# Use awk with a forced tab delimiter (-F'\t') for GFF3 standard compliance.
awk -F'\t' -v chr="$CHR_NAME" -v start="$START_COORD" -v end="$END_COORD" '
# 1. Filter by Chromosome Name
$1 == chr &&
# 2. Filter by feature type "gene"
$3 == "CDS" &&
# 3. CRITICAL: Filter by coordinates (start and end must be INSIDE the inversion boundaries)
# A gene is inside if its END ($5) is before the inversion START ($4 < start)
# OR if its START ($4) is after the inversion END ($4 > end).
( $5 < start || $4 > end ) {

    # 4. Extract and clean the Gene ID from the 9th column (Attributes)
    if (match($9, /ID=([^;]+)/)) {
        # RSTART is the start position of the match, RLENGTH is the length.
        # We start 3 characters after RSTART to skip the "ID=" prefix
        gene_id = substr($9, RSTART + 3, RLENGTH - 3)
        print gene_id
    }
}' "$GFF3_FILE" > "$CDS_LIST"

# --- Report Results ---
COUNT=$(wc -l < "$CDS_LIST")

if [ $? -eq 0 ]; then
    echo "Successfully extracted ${COUNT} gene IDs inside the inversion."
    echo "Gene list saved to: ${CDS_LIST}"
    echo "--- Extraction Complete ---"
else
    echo "ERROR: Awk command failed. Check your GFF3 file format."
fi

####RUN2:ortholog list####


conda activate biopython

python run2_ortholog_mapper.py \
  --gene_list $GENE_LIST \
  --synteny_hits  $SYNTENY_HITS\
  --output_file target_1JMC18_orthologs.txt


####RUN3: extract cds####

python run3_extract_cds.py \
  --genome $GENOME_4JMC19C_FA \
  --gff $GFF3_FILE \
  --genes $GENE_LIST \
  --out 4JMC19C_ortholog_cds.fasta

python run3_extract_cds.py \
  --genome $GENOME_1JMC18_FA \
  --gff $GFF3_FILE_1JMC18 \
  --genes target_1JMC18_orthologs.txt \
  --out 1JMC18_ortholog_cds.fasta

conda deactivate


####RUN4: align orthologs####

conda activate biopython_mafft

python run4_align_orthologs.py \
    4JMC19C_ortholog_cds.fasta \
    1JMC18_ortholog_cds.fasta \
    $SYNTENY_HITS \
    aligned_pairs

####run5: prepare dnds###
python run5_prepare_dnds_data.py \
    aligned_pairs \
    filtered_alignments \
    final_ortholog_list.txt \
    alignment_gap_report.csv


####run6: pal2nal####


python run6_pal2nal.py \
    --pep1 /home/igomez/divergence_rearrangements_helodes/proteins/braker_4JMC19C.aa \
    --pep2 /home/igomez/divergence_rearrangements_helodes/proteins/braker_1JMC18.aa \
    --map $SYNTENY_HITS \
    --orthos final_ortholog_list.txt \
    --filter_dir filtered_alignments \
    --out_dir codon_aligned \
    --pep_tmp temp_pep

python run7_prepare_paml.py \
  --input_dir codon_aligned \
  --output_dir paml_inputs \
  --tree "(Carex_bobo:0.1,Carex_bory:0.1);" \
  --tree_name twoseq.tree

conda deactivate

cd ./paml_inputs
conda activate paml_env

echo "--- Starting PAML Batch Run ---"
echo "Fixing EOF newline issue for all .phy files..."
# Run the newline fix again just to be safe
sed -i '$a\' *.phy

# Loop through all control files (*.ctl) in the current directory
for ctl_file in *.ctl; do
    # Extract the file name (without path) for logging
    base_name=$(basename "$ctl_file")

    # Define a log file for this specific run
    log_file="${base_name%.ctl}.log"

    echo "Processing $base_name -> Output logged to $log_file"

    # Run codeml, redirecting stdout and stderr to the log file.
    # This prevents the terminal mixing that caused the previous error.
    codeml "$ctl_file" > "$log_file" 2>&1

    # Check the exit status of the codeml run
    if [ $? -eq 0 ]; then
        echo "  [SUCCESS] $base_name completed."
    else
        echo "  [FAILURE] $base_name failed. Check $log_file for details."
    fi

done
conda deactivate
