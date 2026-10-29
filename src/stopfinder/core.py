"""Stop-codon analysis of a female/male isoform pair. Used by both transcript and genome mode.

Terms (0-based indices internally, 1-based positions in the output):
  F, M              female and male mRNA, 5'->3' in transcript orientation
  segments          female-specific segments of F: list of (start, end), 0-based, end exclusive
  f2m               for each position of F, the aligned position of M, or None
  female_start      index in F of the first base of the start codon shared with M
  male_start        index in M of the first base of its start codon
"""
from collections import Counter
from dataclasses import dataclass, field

from Bio.Data import CodonTable
from Bio.Seq import Seq

CONTEXT = 6       # nucleotides reported on each side of a stop codon
JUNCTION = 30     # nucleotides on each side of a splice junction (60-nt junction probes)


def stop_codons(table=1):
    return set(CodonTable.unambiguous_dna_by_id[table].stop_codons)


def translate(seq, table=1):
    return str(Seq(seq[: len(seq) // 3 * 3]).translate(table=table))


def rna(codon):
    return codon.replace("T", "U")


def first_stop(seq, start, stops):
    """Index of the first in-frame stop codon at or after `start`, or None."""
    for j in range(start, len(seq) - 2, 3):
        if seq[j:j + 3] in stops:
            return j
    return None


def longest_orf(seq, stops, table=1):
    """Longest ATG-initiated ORF that ends in a stop codon, on the given strand.

    Returns (start, end) with end exclusive and including the stop, or None.
    """
    best = None
    for frame in range(3):
        start = None
        for j in range(frame, len(seq) - 2, 3):
            c = seq[j:j + 3]
            if start is None and c == "ATG":
                start = j
            elif start is not None and c in stops:
                if best is None or j + 3 - start > best[1] - best[0]:
                    best = (start, j + 3)
                start = None
    return best


@dataclass
class PairInput:
    pair_id: str
    female_id: str
    male_id: str
    female: str
    male: str
    segments: list
    f2m: list
    female_start: int
    male_start: int
    mode: str = "transcript"
    female_coord: object = None        # callable: F index -> genomic coordinate string (genome mode)
    segment_coord: object = None       # callable: (start, end) -> genomic interval string (genome mode)
    notes: list = field(default_factory=list)


def analyse(p, table=1):
    """Return (summary_row, per_stop_rows, junction_rows, (female_protein, male_protein)) for one pair."""
    stops = stop_codons(table)
    F, M, fs, ms = p.female, p.male, p.female_start, p.male_start
    notes = list(p.notes)
    spec = [False] * len(F)
    for a, b in p.segments:
        for i in range(a, b):
            spec[i] = True

    fstop = first_stop(F, fs, stops)
    mstop = first_stop(M, ms, stops)
    fprot = translate(F[fs:fstop + 3], table) if fstop is not None else translate(F[fs:], table)
    mprot = translate(M[ms:mstop + 3], table) if mstop is not None else translate(M[ms:], table)
    shared_n = 0
    while shared_n < min(len(fprot), len(mprot)) and fprot[shared_n] == mprot[shared_n]:
        shared_n += 1
    if F[fs:fs + 3] != "ATG":
        notes.append(f"start codon is {F[fs:fs + 3]}, not ATG")
    if fstop is None:
        notes.append("no in-frame stop before the end of the female sequence")
    if mstop is None:
        notes.append("no in-frame stop before the end of the male sequence")

    # the segment the stop belongs to: the last coding-region segment that starts at or before the stop
    # (female-specific 5' UTR sequence, before the start codon, cannot carry the stop)
    cds_segs = [s for s in p.segments if s[1] > fs]
    key = None
    if fstop is not None:
        before = [s for s in cds_segs if s[0] <= fstop]
        key = before[-1] if before else None
    elif cds_segs:
        key = cds_segs[0]

    location = frame = ""
    nt4 = ctx = nxt_n = nxt_c = tail = ""
    if fstop is not None:
        n_spec = sum(spec[fstop:fstop + 3])
        if n_spec == 3:
            location = "female-specific"
        elif n_spec > 0:
            location = "junction"
        elif all(p.f2m[i] is None for i in range(fstop, min(fstop + 3, len(p.f2m)))):
            location = "not aligned to male (UTR, contig end or small indel)"
        elif not cds_segs:
            location = "shared (no female-specific sequence)"
        elif fstop < cds_segs[0][0]:
            location = "shared, upstream of female-specific sequence"
            notes.append("translation ends before the female-specific sequence, which is therefore not coding")
        else:
            location = "shared, downstream of female-specific sequence"
        pm = p.f2m[fstop] if fstop < len(p.f2m) else None
        if pm is not None:
            frame = "male frame" if (pm - ms) % 3 == 0 else f"shifted +{(pm - ms) % 3}"
        if location == "shared, downstream of female-specific sequence" and key is not None:
            seg_len = key[1] - key[0]
            if frame.startswith("shifted"):
                notes.append(f"no in-frame stop in the female-specific segment ({seg_len} bp, "
                             f"{'not ' if seg_len % 3 else ''}a multiple of 3); the reading frame is shifted "
                             f"and the stop lies {fstop - key[1]} nt into the downstream shared sequence")
            elif frame == "male frame":
                notes.append(f"no in-frame stop in the female-specific segment ({seg_len} bp); translation "
                             f"continues in the male frame")
        nt4 = F[fstop + 3] if fstop + 3 < len(F) else ""
        ctx = F[max(fstop - CONTEXT, 0):fstop + 3 + CONTEXT]
        nxt = first_stop(F, fstop + 3, stops)
        if nxt is not None:
            nxt_n, nxt_c = (nxt - fstop) // 3, F[nxt:nxt + 3]
        if key is not None:
            first = fs + max(key[0] - fs, 0) // 3 * 3
            tail = translate(F[first:fstop + 3], table)

    # every stop codon inside the female-specific segments, in all three frames
    per_stop, inframe = [], Counter()
    for si, (a, b) in enumerate(p.segments, 1):
        for j in range(a, b - 2):
            c = F[j:j + 3]
            if c in stops:
                fr = (j - fs) % 3
                if fr == 0 and j >= fs:
                    inframe[c] += 1
                per_stop.append({
                    "pair_id": p.pair_id, "segment": si, "frame_vs_orf": fr,
                    "codon_index": (j - fs) // 3 + 1 if fr == 0 and j >= fs else "",
                    "codon": c, "codon_rna": rna(c), "pos_female": j + 1,
                    "genomic": p.female_coord(j) if p.female_coord else "",
                    "is_termination_codon": j == fstop,
                    "context": F[max(j - CONTEXT, 0):j + 3 + CONTEXT],
                })

    def seg_txt(a, b, mark):
        g = f" {p.segment_coord(a, b)}" if p.segment_coord else ""
        return f"{a + 1}-{b}({b - a}bp){g}{'*' if mark else ''}"

    # splice-junction probes for the segment that carries the stop
    junctions = []
    if key is not None:
        a = key[0]
        if a >= JUNCTION and a + JUNCTION <= len(F):
            junctions.append({"pair_id": p.pair_id, "junction": "female", "sequence": F[a - JUNCTION:a + JUNCTION]})
        pm = p.f2m[a - 1] if a > 0 else None
        if pm is not None and pm + 1 >= JUNCTION and pm + 1 + JUNCTION <= len(M):
            junctions.append({"pair_id": p.pair_id, "junction": "male", "sequence": M[pm + 1 - JUNCTION:pm + 1 + JUNCTION]})

    row = {
        "pair_id": p.pair_id, "mode": p.mode, "female_id": p.female_id, "male_id": p.male_id,
        "female_cds": f"{fs + 1}-{fstop + 3 if fstop is not None else ''}",
        "male_cds": f"{ms + 1}-{mstop + 3 if mstop is not None else ''}",
        "start_codon": F[fs:fs + 3],
        "female_protein_aa": len(fprot.rstrip("*")), "male_protein_aa": len(mprot.rstrip("*")),
        "shared_nterm_aa": shared_n,
        "female_specific_segments": ";".join(seg_txt(a, b, (a, b) == key) for a, b in p.segments),
        "female_specific_bp": key[1] - key[0] if key else "",
        "stop_codon": F[fstop:fstop + 3] if fstop is not None else "",
        "stop_codon_rna": rna(F[fstop:fstop + 3]) if fstop is not None else "",
        "stop_pos_female": fstop + 1 if fstop is not None else "",
        "stop_genomic": p.female_coord(fstop) if p.female_coord and fstop is not None else "",
        "stop_location": location, "frame_vs_male": frame,
        "nt_plus4": nt4, "context": ctx,
        "codons_to_next_stop": nxt_n, "next_stop_codon": nxt_c, "next_stop_codon_rna": rna(nxt_c),
        "inframe_stops_in_female_specific": ",".join(f"{k}:{v}" for k, v in sorted(inframe.items())) or "none",
        "female_tail_peptide": tail,
        "note": "; ".join(notes),
    }
    proteins = (fprot, mprot)
    return row, per_stop, junctions, proteins


SUMMARY_COLUMNS = [
    "pair_id", "mode", "female_id", "male_id", "female_cds", "male_cds", "start_codon",
    "female_protein_aa", "male_protein_aa", "shared_nterm_aa", "female_specific_segments", "female_specific_bp",
    "stop_codon", "stop_codon_rna", "stop_pos_female", "stop_genomic", "stop_location", "frame_vs_male",
    "nt_plus4", "context", "codons_to_next_stop", "next_stop_codon", "next_stop_codon_rna",
    "inframe_stops_in_female_specific", "female_tail_peptide", "note",
]
STOP_COLUMNS = ["pair_id", "segment", "frame_vs_orf", "codon_index", "codon", "codon_rna", "pos_female",
                "genomic", "is_termination_codon", "context"]
JUNCTION_COLUMNS = ["pair_id", "junction", "sequence"]
