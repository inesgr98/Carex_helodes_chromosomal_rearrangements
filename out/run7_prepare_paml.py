#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare PHYLIP + ctl files for PAML")

    parser.add_argument("--input_dir", required=True,
                        help="Directory containing codon-aligned FASTA (.codon.fasta)")
    parser.add_argument("--output_dir", required=True,
                        help="Where PHYLIP and ctl files will be written")
    parser.add_argument("--tree", required=True,
                        help="Tree string for PAML (e.g. '(A:0.1,B:0.1);')")
    parser.add_argument("--tree_name", default="twoseq.tree",
                        help="Filename for the tree file")

    return parser.parse_args()

def write_paml_phylip(records, output_path):
    """Writes PAML-compatible PHYLIP sequential file."""
    if not records:
        print(f"Warning: No records for {output_path}")
        return

    num_seqs = len(records)
    aln_length = len(records[0].seq)
    block_size = 10

    with open(output_path, 'w') as f:
        f.write(f" {num_seqs} {aln_length}\n")
        for rec in records:
            formatted_name = f"{rec.id:<10}"
            f.write(f"{formatted_name}  ")

            for start in range(0, aln_length, block_size):
                block = str(rec.seq[start:start + block_size])
                if start == 0:
                    f.write(block)
                else:
                    f.write(f" {block}")
            f.write("\n")

def generate_control_file(ctl_path, aln_path, tree_path, gene_id):
    """Writes a basic codeml control file (Model 0)."""
    ctl_content = f"""
seqfile = {os.path.basename(aln_path)}
outfile = {gene_id}.mlc
treefile = {os.path.basename(tree_path)}

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

def prepare_paml_input(args):
    """Main pipeline."""
    os.makedirs(args.output_dir, exist_ok=True)

    # Write tree file
    tree_path = os.path.join(args.output_dir, args.tree_name)
    with open(tree_path, 'w') as f:
        f.write(args.tree)
    print(f"Created tree file: {tree_path}")

    success_count = 0

    for filename in os.listdir(args.input_dir):
        if filename.endswith(".codon.fasta"):
            gene_id = filename.replace(".codon.fasta", "")
            fasta_path = os.path.join(args.input_dir, filename)
            phylip_path = os.path.join(args.output_dir, f"{gene_id}.phy")
            ctl_path = os.path.join(args.output_dir, f"{gene_id}.ctl")

            try:
                records = list(SeqIO.parse(fasta_path, "fasta"))

                if len(records) != 2:
                    print(f"ERROR: {filename} has {len(records)} seqs; expected 2. Skipping.")
                    continue

                # Force names to match tree tips
                records[0].id = records[0].name = "Carex_bobo"
                records[0].description = ""
                records[1].id = records[1].name = "Carex_bory"
                records[1].description = ""

                write_paml_phylip(records, phylip_path)
                generate_control_file(ctl_path, phylip_path, tree_path, gene_id)

                success_count += 1

            except Exception as e:
                print(f"FAILED {gene_id}: {e}")

    print("-" * 40)
    print(f"Prepared {success_count} PAML input files in '{args.output_dir}/'.")

if __name__ == "__main__":
    args = parse_args()
    prepare_paml_input(args)

