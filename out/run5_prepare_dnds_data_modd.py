#!/usr/bin/env python3
import os
from Bio import SeqIO
import csv

# --- Configuration ---
INPUT_DIR = "aligned_pairs"
OUTPUT_DIR = "filtered_alignments"
FINAL_LIST_FILE = "final_ortholog_list_chr35.txt"
GAP_REPORT_FILE = "alignment_gap_report_chr35.csv"

MIN_ALIGNMENT_LENGTH_BP = 120   # Minimum nucleotide length (30 codons)
MAX_GAP_PERCENTAGE = 0.10      # Maximum allowed gap percentage (10%)
# ---------------------------------------------------------

def calculate_gap_percentage(sequence):
    """Calculates the fraction of gaps ('-' or '?') in a sequence."""
    total_length = len(sequence)
    if total_length == 0:
        return 1.0
    gaps = sequence.count('-') + sequence.count('?')
    return gaps / total_length


def filter_and_prepare_alignments():
    """Filters alignments by length and gap content, and writes passing ones."""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    all_files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".aln.fasta")]
    total = len(all_files)
    passed = 0
    passed_ids = []
    gap_summary = []  # for CSV output

    print(f"Checking {total} alignments from '{INPUT_DIR}/'...")

    for filename in sorted(all_files):
        input_path = os.path.join(INPUT_DIR, filename)
        cds_id = filename.replace(".aln.fasta", "")

        try:
            records = list(SeqIO.parse(input_path, "fasta"))
            if len(records) != 2:
                print(f"Skipping {cds_id}: not exactly 2 sequences.")
                continue

            seq1, seq2 = str(records[0].seq).upper(), str(records[1].seq).upper()

            # --- 1️⃣ Length filter ---
            if len(seq1) < MIN_ALIGNMENT_LENGTH_BP:
                print(f"Skipping {cds_id}: too short ({len(seq1)} bp).")
                continue

            # --- 2️⃣ Gap detection ---
            gap1 = calculate_gap_percentage(seq1)
            gap2 = calculate_gap_percentage(seq2)
            avg_gap = (gap1 + gap2) / 2

            has_gaps = ("YES" if ("-" in seq1 or "-" in seq2) else "NO")

            if gap1 > MAX_GAP_PERCENTAGE or gap2 > MAX_GAP_PERCENTAGE:
                passed_flag = "NO"
                print(f"Skipping {cds_id}: too many gaps (Seq1={gap1:.2%}, Seq2={gap2:.2%}).")
            else:
                passed_flag = "YES"
                # --- 3️⃣ Write passing alignment ---
                output_path = os.path.join(OUTPUT_DIR, filename)
                SeqIO.write(records, output_path, "fasta")
                passed_ids.append(cds_id)
                passed += 1

            # Store for gap report
            gap_summary.append({
                "CDS_ID": cds_id,
                "Seq1_gap_%": f"{gap1*100:.2f}",
                "Seq2_gap_%": f"{gap2*100:.2f}",
                "Has_gaps": has_gaps,
                "Avg_gap_%": f"{avg_gap*100:.2f}",
                "Passed_filter": passed_flag
            })

        except Exception as e:
            print(f"Error processing {cds_id}: {e}")

    # --- Write list of passing CDS ---
    with open(FINAL_LIST_FILE, "w") as f:
        f.write("\n".join(passed_ids))

    # --- Write CSV gap report ---
    with open(GAP_REPORT_FILE, "w", newline="") as csvfile:
        fieldnames = ["CDS_ID", "Seq1_gap_%", "Seq2_gap_%", "Avg_gap_%", "Has_gaps", "Passed_filter"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(gap_summary)

    print("\n--- Filtering Summary ---")
    print(f"Total input alignments: {total}")
    print(f"Alignments passed: {passed}")
    print(f"Filtered files saved to: {OUTPUT_DIR}/")
    print(f"List of passing CDS saved to: {FINAL_LIST_FILE}")
    print(f"Gap report saved to: {GAP_REPORT_FILE}")
    print("-" * 40)


if __name__ == "__main__":
    print("--- Starting CDS Alignment Filtering for dN/dS ---")
    filter_and_prepare_alignments()
    print("Next step: PAL2NAL or codon-based dN/dS analysis.")

