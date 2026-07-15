import os
from Bio import SeqIO

# --- Configuration ---
INPUT_DIR = "codon_aligned"   # CODON-ALIGNED CDS OUTPUT
OUTPUT_DIR = "paml_inputs"    # Where PHYLIP + ctl files will go

# Standard PAML tree for a two-sequence comparison (Borbonica and Boryana)
TREE_FILE_CONTENT = "(Carex_bobo:0.1, Carex_bory:0.1);"
TREE_FILE_NAME = "twoseq.tree"

# ---------------------

def write_paml_phylip(records, output_path):
    """
    Writes a PAML-compatible PHYLIP sequential file.
    Names are 10 characters max; sequences are broken into blocks of 10.
    """
    if not records:
        print(f"Warning: No records to write for {output_path}")
        return

    num_seqs = len(records)
    aln_length = len(records[0].seq)
    block_size = 10

    with open(output_path, 'w') as f:
        f.write(f" {num_seqs} {aln_length}\n")
        for rec in records:
            formatted_name = f"{rec.id:<10}"
            seq_block = str(rec.seq[:block_size])
            f.write(f"{formatted_name}  {seq_block}")
            for start in range(block_size, aln_length, block_size):
                seq_block = str(rec.seq[start:start + block_size])
                f.write(f" {seq_block}")
            f.write("\n")

def generate_control_file(ctl_path, aln_path, tree_path, gene_id):
    """Generates a simple codeml control file (Model 0) for pairwise dN/dS."""
    ctl_content = f"""
seqfile = {os.path.basename(aln_path)}
outfile = {gene_id}.mlc
treefile = {os.path.basename(tree_path)}

# Model parameters
noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0

model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2.0
fix_omega = 0
omega = 0.5

Small_Diff = 0.5e-6
cleandata = 1
"""
    with open(ctl_path, 'w') as f:
        f.write(ctl_content.strip())

def prepare_paml_input():
    """Main pipeline for converting codon-aligned FASTA to PHYLIP + ctl for PAML."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Write shared tree
    tree_path = os.path.join(OUTPUT_DIR, TREE_FILE_NAME)
    with open(tree_path, 'w') as f:
        f.write(TREE_FILE_CONTENT)
    print(f"Created shared tree file: {tree_path}")

    success_count = 0

    for filename in os.listdir(INPUT_DIR):
        if filename.endswith(".codon.fasta"):
            gene_id = filename.replace(".codon.fasta", "")
            fasta_input_path = os.path.join(INPUT_DIR, filename)
            phylip_output_path = os.path.join(OUTPUT_DIR, f"{gene_id}.phy")
            ctl_output_path = os.path.join(OUTPUT_DIR, f"{gene_id}.ctl")

            try:
                records = list(SeqIO.parse(fasta_input_path, "fasta"))

                if len(records) != 2:
                    print(f"ERROR: {filename} does not contain exactly two sequences ({len(records)} found). Skipping.")
                    continue

                # Force 10-character names to match the tree
                records[0].id = records[0].name = "Carex_bobo"
                records[0].description = ""
                records[1].id = records[1].name = "Carex_bory"
                records[1].description = ""

                write_paml_phylip(records, phylip_output_path)
                generate_control_file(ctl_output_path, phylip_output_path, tree_path, gene_id)
                success_count += 1

            except Exception as e:
                print(f"FAILED to process {gene_id}: {e}")

    print("-" * 40)
    print(f"Successfully prepared {success_count} gene inputs in '{OUTPUT_DIR}/'")

if __name__ == "__main__":
    prepare_paml_input()

