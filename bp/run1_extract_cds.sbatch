#!/bin/bash
#SBATCH --job-name=extract_genes
#SBATCH --output=extract_genes_%j.log
#SBATCH --error=extract_genes_%j.err
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
# --- Configuration ---
GFF3_FILE="/home/igomez/genome_assembly_helodes_thesis/C2b_braker/4JMC19C/braker/braker.gff3"
OUTPUT_FILE="target_cds_chr4.txt" # Updated output file name
CHR_NAME="scaffold_16"

# Inversion coordinates for Carex borbonica
START_COORD="4878165"
END_COORD="9622376"

echo "--- Starting Gene ID Extraction (Outside Inversion) ---"
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
( $4 >= start && $5 <= end ) {

    # 4. Extract and clean the Gene ID from the 9th column (Attributes)
    if (match($9, /ID=([^;]+)/)) {
        # RSTART is the start position of the match, RLENGTH is the length.
        # We start 3 characters after RSTART to skip the "ID=" prefix
        gene_id = substr($9, RSTART + 3, RLENGTH - 3)
        print gene_id
    }
}' "$GFF3_FILE" > "$OUTPUT_FILE"

# --- Report Results ---
COUNT=$(wc -l < "$OUTPUT_FILE")

if [ $? -eq 0 ]; then
    echo "Successfully extracted ${COUNT} gene IDs inside the inversion."
    echo "Gene list saved to: ${OUTPUT_FILE}"
    echo "--- Extraction Complete ---"
else
    echo "ERROR: Awk command failed. Check your GFF3 file format."
fi

echo "---"
echo "If the count is 0, the scaffold name might be wrong, or your GFF3 is incorrectly formatted."
echo "NEXT STEP: Run 'run2_ortholog_mapper.py' (Step 2) using this new file."

