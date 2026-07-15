#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
import csv


def calculate_gap_percentage(sequence):
    """Calculates fraction of gaps ('-' or '?') in a sequence."""
    total_length = len(sequence)
    if total_length == 0:
        return 1.0
    gaps = sequence.count('-') + sequence.count('?')
    return gaps / total_length


def filter_and_prepare_alignments(input_dir, output_dir, final_list_file, gap_report_file,
                                  min_alignment_length_bp=120, max_gap_percentage=0.10):
    """Filters alignments by length and gap content, and writes passing ones."""

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_files = [f for f in os.listdir(input_dir) if f.endswith(".aln.fasta")]
    total = len(all_files)
    passed = 0
    passed_ids = []
    gap_summary = []

    print(f"Checking {total} alignments from '{input_dir}/'...")

    for filename in sorted(all_files):
        input_path = os.path.join(input_dir, filename)
        cds_id = filename.replace(".aln.fasta", "")

        try:
            records = list(SeqIO.parse(input_path, "fasta"))
            if len(records) != 2:
                print(f"Skipping {cds_id}: not exactly 2 sequences.")
                continue

            seq1 = str(records[0].seq).upper()
            seq2 = str(records[1].seq).upper()

            # ----- 1️⃣ Length filter -----
            if len(seq1) < min_alignment_length_bp:
                print(f"Skipping {cds_id}: too short ({len(seq1)} bp).")
                continue

            # ----- 2️⃣ Gap detection -----
            gap1 = calculate_gap_percentage(seq1)
            gap2 = calculate_gap_percentage(seq2)
            avg_gap = (gap1 + gap2) / 2

            has_gaps = ("YES" if ("-" in seq1 or "-" in seq2) else "NO")

            if gap1 > max_gap_percentage or gap2 > max_gap_percentage:
                passed_flag = "NO"
                print(f"Skipping {cds_id}: too many gaps (Seq1={gap1:.2%}, Seq2={gap2:.2%}).")
            else:
                passed_flag = "YES"
                output_path = os.path.join(output_dir, filename)
                SeqIO.write(records, output_path, "fasta")
                passed_ids.append(cds_id)
                passed += 1

            # Add to gap report
            gap_summary.append({
                "CDS_ID": cds_id,
                "Seq1_gap_%": f"{gap1*100:.2f}",
                "Seq2_gap_%": f"{gap2*100:.2f}",
                "Avg_gap_%": f"{avg_gap*100:.2f}",
                "Has_gaps": has_gaps,
                "Passed_filter": passed_flag
            })

        except Exception as e:
            print(f"Error processing {cds_id}: {e}")

    # ----- Save passing CDS list -----
    with open(final_list_file, "w") as f:
        f.write("\n".join(passed_ids))

    # ----- Save CSV gap report -----
    with open(gap_report_file, "w", newline="") as csvfile:
        fieldnames = ["CDS_ID", "Seq1_gap_%", "Seq2_gap_%", "Avg_gap_%",
                      "Has_gaps", "Passed_filter"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(gap_summary)

    print("\n--- Filtering Summary ---")
    print(f"Total input alignments: {total}")
    print(f"Alignments passed: {passed}")
    print(f"Filtered files saved to: {output_dir}/")
    print(f"List of passing CDS saved to: {final_list_file}")
    print(f"Gap report saved to: {gap_report_file}")
    print("-" * 40)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("\nUsage:")
        print("  python run5_filter_alignments.py INPUT_DIR OUTPUT_DIR final_list.txt gap_report.csv\n")
        sys.exit(1)

    INPUT_DIR = sys.argv[1]
    OUTPUT_DIR = sys.argv[2]
    FINAL_LIST_FILE = sys.argv[3]
    GAP_REPORT_FILE = sys.argv[4]

    print("--- Starting RUN5: CDS Alignment Filtering ---")

    filter_and_prepare_alignments(
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        final_list_file=FINAL_LIST_FILE,
        gap_report_file=GAP_REPORT_FILE
    )

    print("RUN5 complete. Next step: PAL2NAL / codon alignment for dN/dS.")

