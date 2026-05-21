#!/usr/bin/env python3
import os
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(description="Run PAL2NAL codon alignment workflow")

    parser.add_argument("--pep1", required=True,
                        help="Peptide FASTA for species 1")
    parser.add_argument("--pep2", required=True,
                        help="Peptide FASTA for species 2")
    parser.add_argument("--map", required=True,
                        help="Ortholog mapping file (.synHits.txt)")
    parser.add_argument("--orthos", required=True,
                        help="Final ortholog list (CDS IDs)")
    parser.add_argument("--filter_dir", required=True,
                        help="Directory containing filtered nucleotide alignments")
    parser.add_argument("--out_dir", default="codon_aligned",
                        help="Output directory for codon alignments")
    parser.add_argument("--pep_tmp", default="temp_pep",
                        help="Temporary peptide alignment dir")

    return parser.parse_args()

def get_base_gene_id(full_id):
    return full_id.split('|')[0].split('.')[0].strip()

def load_fasta_to_dict(filepath):
    fasta_dict = {}
    for rec in SeqIO.parse(filepath, "fasta"):
        base_id = get_base_gene_id(rec.id)
        fasta_dict[base_id] = rec
    return fasta_dict

def load_ortholog_map(map_file):
    mapping = {}
    with open(map_file, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith("ofID1"):
                continue
            parts = line.strip().split()
            if len(parts) >= 13:
                id_borbonica = get_base_gene_id(parts[12])
                id_boryana = get_base_gene_id(parts[4])
                mapping[id_borbonica] = id_boryana
    return mapping

def run_pal2nal_for_cds(args):
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.pep_tmp, exist_ok=True)

    pep_dict_1 = load_fasta_to_dict(args.pep1)
    pep_dict_2 = load_fasta_to_dict(args.pep2)
    ortholog_map = load_ortholog_map(args.map)

    with open(args.orthos, 'r') as f:
        cds_ids = [line.strip() for line in f if line.strip()]

    total = len(cds_ids)
    print(f"Processing {total} CDS alignments for PAL2NAL...")

    for i, cds_id in enumerate(cds_ids, start=1):
        base_gene = cds_id.split('_CDS')[0]
        partner_gene = ortholog_map.get(base_gene)

        if not partner_gene:
            print(f"[{i}/{total}] Skipping {cds_id}: partner not found.")
            continue

        nuc_align_path = os.path.join(args.filter_dir, f"{cds_id}.aln.fasta")
        if not os.path.exists(nuc_align_path):
            print(f"[{i}/{total}] Missing alignment: {cds_id}")
            continue

        nucleotide_records = list(SeqIO.parse(nuc_align_path, "fasta"))
        protein_records = []

        for rec in nucleotide_records:
            clean_seq = Seq(str(rec.seq).replace('-', '').replace('?', ''))
            try:
                prot_seq = clean_seq.translate(to_stop=True)
                protein_records.append(SeqRecord(prot_seq, id=rec.id, description=""))
            except Exception as e:
                print(f"Translation error for {rec.id}: {e}")
                continue

        pep_out_path = os.path.join(args.pep_tmp, f"{cds_id}.pep.fasta")
        SeqIO.write(protein_records, pep_out_path, "fasta")

        pal2nal_out_path = os.path.join(args.out_dir, f"{cds_id}.codon.fasta")
        pal2nal_cmd = [
            "pal2nal.pl",
            pep_out_path,
            nuc_align_path,
            "-output", "fasta",
            "-nogap"
        ]

        try:
            result = subprocess.run(pal2nal_cmd, check=True, capture_output=True, text=True)
            with open(pal2nal_out_path, 'w') as f:
                f.write(result.stdout)
            print(f"[{i}/{total}] SUCCESS: {cds_id}")
        except subprocess.CalledProcessError as e:
            print(f"[{i}/{total}] PAL2NAL failed: {cds_id}: {e.stderr}")

if __name__ == "__main__":
    args = parse_args()
    print("--- Starting PAL2NAL codon alignment ---")
    run_pal2nal_for_cds(args)

