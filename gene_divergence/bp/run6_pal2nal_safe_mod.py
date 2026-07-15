#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# --- Configuration ---
PEP_FASTA_1 = "/home/igomez/divergence_rearrangements_helodes/proteins/braker_4JMC19C.aa"
PEP_FASTA_2 = "/home/igomez/divergence_rearrangements_helodes/proteins/braker_1JMC18.aa"
MAP_FILE = "/home/igomez/genome_assembly_helodes_thesis/C3b_syntheny/syntenicHits/C_helodes_4JMC19C_vs_C_helodes_1JMC18.synHits.txt"

FILTERED_NUCLEOTIDE_DIR = "filtered_alignments"
FINAL_ORTHOLOG_LIST = "final_ortholog_list_chr35.txt"

# New output directory for codon-aligned files
OUTPUT_DIR = "codon_aligned"
TEMP_PEP_DIR = "temp_pep"
# ----------------------

def get_base_gene_id(full_id):
    """Extracts the base ID (e.g. g28319) from a FASTA header."""
    return full_id.split('|')[0].split('.')[0].strip()

def load_fasta_to_dict(filepath):
    """Loads a FASTA file as { base_gene_id : SeqRecord }."""
    fasta_dict = {}
    for rec in SeqIO.parse(filepath, "fasta"):
        base_id = get_base_gene_id(rec.id)
        fasta_dict[base_id] = rec
    return fasta_dict

def load_ortholog_map(map_file):
    """Loads ortholog map as { borbonica_ID : boryana_ID }."""
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

def run_pal2nal_for_cds():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(TEMP_PEP_DIR, exist_ok=True)

    pep_dict_1 = load_fasta_to_dict(PEP_FASTA_1)
    pep_dict_2 = load_fasta_to_dict(PEP_FASTA_2)
    ortholog_map = load_ortholog_map(MAP_FILE)

    with open(FINAL_ORTHOLOG_LIST, 'r') as f:
        cds_ids = [line.strip() for line in f if line.strip()]  # e.g. g28319_CDS1

    total = len(cds_ids)
    print(f"Processing {total} CDS alignments for PAL2NAL...")

    for i, cds_id in enumerate(cds_ids, start=1):
        # Get gene base (e.g. g28319 from g28319_CDS1)
        base_gene = cds_id.split('_CDS')[0]
        partner_gene = ortholog_map.get(base_gene)

        if not partner_gene:
            print(f"[{i}/{total}] Skipping {cds_id}: partner not found in map.")
            continue

        nuc_align_path = os.path.join(FILTERED_NUCLEOTIDE_DIR, f"{cds_id}.aln.fasta")
        if not os.path.exists(nuc_align_path):
            print(f"[{i}/{total}] Skipping {cds_id}: alignment not found.")
            continue

        # Load aligned nucleotide sequences and translate to protein
        nucleotide_records = list(SeqIO.parse(nuc_align_path, "fasta"))
        protein_records = []
        for rec in nucleotide_records:
            clean_seq = Seq(str(rec.seq).replace('-', '').replace('?', ''))
            try:
                prot_seq = clean_seq.translate(to_stop=True)
                protein_records.append(SeqRecord(prot_seq, id=rec.id, description=""))
            except Exception as e:
                print(f"⚠️  Skipping {rec.id} due to translation error: {e}")
                continue

        pep_out_path = os.path.join(TEMP_PEP_DIR, f"{cds_id}.pep.fasta")
        SeqIO.write(protein_records, pep_out_path, "fasta")

        # Run PAL2NAL
        pal2nal_out_path = os.path.join(OUTPUT_DIR, f"{cds_id}.codon.fasta")
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
            print(f"[{i}/{total}] ✅ SUCCESS: {cds_id}")
        except subprocess.CalledProcessError as e:
            print(f"[{i}/{total}] ❌ PAL2NAL failed for {cds_id}: {e.stderr.strip()}")

    print(f"\n--- Finished ---")
    print(f"Codon-aligned CDS files saved to '{OUTPUT_DIR}/'")
    print(f"Temporary peptide files saved to '{TEMP_PEP_DIR}/'")
    print("You can now use these for codon-based dN/dS analysis (e.g., with PAML or Ka/Ks tools).")

if __name__ == "__main__":
    print("--- Starting PAL2NAL codon alignment for CDS ---")
    run_pal2nal_for_cds()

