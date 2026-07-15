#!/usr/bin/env python3
import os
import re
import subprocess
import sys
from collections import defaultdict


# ------------------------------
#   Helper Functions
# ------------------------------

def read_fasta_by_id(fasta_path):
    sequences = {}
    with open(fasta_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            if line.startswith('>'):
                if header:
                    key = make_cds_key(header)
                    if key:
                        sequences[key] = {"header": header, "seq": "".join(seq_lines)}
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip().upper())

        if header:
            key = make_cds_key(header)
            if key:
                sequences[key] = {"header": header, "seq": "".join(seq_lines)}

    print(f"Loaded {len(sequences)} sequences from {os.path.basename(fasta_path)}.")
    return sequences


def make_cds_key(header_line):
    parts = header_line.replace(">", "").split("|")
    if len(parts) < 3:
        return None
    gene_id = parts[0]
    cds_part = parts[2]
    return f"{gene_id}_{cds_part}"


def normalize_gene_id(gene_id):
    normalized = gene_id.split('|')[0]
    if '.t' in normalized:
        normalized = normalized.split('.t')[0]
    return normalized.strip()


def read_synteny_hits_map(map_path):
    print(f"Reading ortholog map from {map_path}...")
    pairs = {}

    try:
        with open(map_path, 'r') as f:
            header = f.readline()
            for line in f:
                if not line.strip():
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 13:
                    continue

                boryana_raw = cols[12].strip()
                borbonica_raw = cols[4].strip()

                borbonica = normalize_gene_id(borbonica_raw)
                boryana = normalize_gene_id(boryana_raw)

                pairs[borbonica] = boryana

    except FileNotFoundError:
        print(f"FATAL ERROR: Map file not found: {map_path}")
        return None

    print(f"Loaded {len(pairs)} ortholog pairs.")
    return pairs


def run_mafft_alignment(input_path, output_path):
    try:
        subprocess.run(
            ["mafft", "--quiet", "--auto", input_path],
            check=True,
            stdout=open(output_path, "w"),
            stderr=subprocess.PIPE
        )
        return True
    except:
        return False


def align_cds_by_orthologs(borbonica_fasta, boryana_fasta, output_dir, map_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ortholog_map = read_synteny_hits_map(map_file)
    if not ortholog_map:
        return

    print("Reading borbonica CDS...")
    borbonica_seqs = read_fasta_by_id(borbonica_fasta)

    print("Reading boryana CDS...")
    boryana_seqs = read_fasta_by_id(boryana_fasta)

    paired_count = 0
    aligned_count = 0

    print("\nStarting pairwise CDS alignments...")

    for borbonica_gene, boryana_gene in ortholog_map.items():

        cds_bor = [k for k in borbonica_seqs if k.startswith(f"{borbonica_gene}_CDS")]

        for cds_key_bor in cds_bor:
            cds_number = cds_key_bor.split("_CDS")[-1]
            cds_key_bor_other = f"{boryana_gene}_CDS{cds_number}"

            if cds_key_bor_other not in boryana_seqs:
                continue

            paired_count += 1

            bor = borbonica_seqs[cds_key_bor]
            bry = boryana_seqs[cds_key_bor_other]

            temp_fa = os.path.join(output_dir, f"{cds_key_bor}.tmp.fa")
            aln_out = os.path.join(output_dir, f"{cds_key_bor}.aln.fasta")

            with open(temp_fa, "w") as tmp:
                tmp.write(f"{bor['header']} (Borbonica)\n{bor['seq']}\n")
                tmp.write(f"{bry['header']} (Boryana)\n{bry['seq']}\n")

            if run_mafft_alignment(temp_fa, aln_out):
                aligned_count += 1

            os.remove(temp_fa)

    print("\n--- Alignment Summary ---")
    print(f"CDS pairs aligned: {aligned_count} / {paired_count}")
    print(f"Alignments saved in: {output_dir}")
    print("-" * 40)


# ------------------------------
#           MAIN
# ------------------------------
if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("\nUsage:")
        print("  python run4_align.py BORBONICA_CDS.fasta BORYANA_CDS.fasta ORTHOLOG_MAP.txt output_dir\n")
        sys.exit(1)

    BORBONICA_FASTA = sys.argv[1]
    BORYANA_FASTA = sys.argv[2]
    MAP_FILE = sys.argv[3]
    OUTPUT_DIR = sys.argv[4]

    print("\n--- Running RUN3: CDS Alignment ---\n")

    align_cds_by_orthologs(
        borbonica_fasta=BORBONICA_FASTA,
        boryana_fasta=BORYANA_FASTA,
        output_dir=OUTPUT_DIR,
        map_file=MAP_FILE
    )

    print("RUN4 completed successfully.\n")

