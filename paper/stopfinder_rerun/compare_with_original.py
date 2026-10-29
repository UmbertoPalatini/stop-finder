#!/usr/bin/env python
"""Compare stopfinder output with the original tables of the paper (steps 02 and 05).

Usage: python compare_with_original.py ORIGINAL_RESULTS_DIR STOPFINDER_GENOME_OUT STOPFINDER_TRANSCRIPTS_OUT
ORIGINAL_RESULTS_DIR holds 02_ref_stops/ref_pairs.tsv and 05_toxwye/toxwye_pairs.tsv.
Prints one line per pair and field that differs, and exits with status 1 if anything differs.
"""
import csv
import sys

orig_dir, gen_dir, tx_dir = sys.argv[1:4]


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


# original column -> stopfinder column
FIELDS = [("female_stop_codon", "stop_codon"), ("nt4_after_stop", "nt_plus4"),
          ("stop_context_-6_+6", "context"), ("codons_to_next_inframe_stop", "codons_to_next_stop"),
          ("next_stop_codon", "next_stop_codon"), ("female_prot_len", "female_protein_aa"),
          ("male_prot_len", "male_protein_aa"), ("female_specific_inframe_stops", "inframe_stops_in_female_specific")]

old = {}
for r in rows(f"{orig_dir}/02_ref_stops/ref_pairs.tsv"):
    old[f"{r['species']}_{r['gene']}"] = r
for r in rows(f"{orig_dir}/05_toxwye/toxwye_pairs.tsv"):
    old[f"{r['species']}_{r['gene']}_{r['isoform']}"] = r
new = {r["pair_id"]: r for d in (gen_dir, tx_dir) for r in rows(f"{d}/summary.tsv")}

diffs = 0
for pid in sorted(old):
    if pid not in new:
        print(f"{pid}\tmissing from stopfinder output")
        diffs += 1
        continue
    o, n = old[pid], new[pid]
    for oc, nc in FIELDS:
        ov, nv = o.get(oc, ""), n.get(nc, "")
        if oc == "female_specific_inframe_stops":
            ov = ov or "none"
        if ov != nv:
            print(f"{pid}\t{nc}\toriginal={ov}\tstopfinder={nv}")
            diffs += 1
    o_in = o["stop_in_female_specific"] == "True"
    n_in = n["stop_location"] == "female-specific"
    if o_in != n_in:
        print(f"{pid}\tstop_in_female_specific\toriginal={o_in}\tstopfinder={n['stop_location']}")
        diffs += 1
print(f"{len(old)} pairs compared, {diffs} difference(s)")
sys.exit(1 if diffs else 0)
