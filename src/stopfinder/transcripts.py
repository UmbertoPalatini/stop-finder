"""Transcript mode: a female and a male mRNA (any source) are aligned to each other.

The female-specific sequence is the part of the female mRNA that does not align to the male mRNA.
"""
import re

from Bio.Align import PairwiseAligner

from .core import PairInput, longest_orf, stop_codons

_CDS_TAG = re.compile(r"\bCDS[=:](\d+)-(\d+)")


def make_aligner():
    # Global alignment with free end gaps, scored for near-identical sequence (two isoforms of one gene):
    # a mismatch costs 4 matches and long gaps are cheap to extend. Shared exons align despite a few SNPs
    # between individuals, while an inserted or replaced exon stays unaligned instead of being forced
    # into a poor alignment.
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 1, -4
    a.open_gap_score, a.extend_gap_score = -10, -0.1
    a.end_gap_score = 0
    return a


def align(F, M, min_gap=15, min_island=30):
    """Return f2m (female index -> male index or None) and the identity of the aligned columns.

    Aligned blocks separated by small gaps (< min_gap in both sequences: indels between individuals) are
    grouped into islands. Islands shorter than min_island nt are chance matches between unrelated sequence
    (for example two different terminal exons) and are discarded.
    """
    aln = make_aligner().align(F, M)[0]
    blocks = [((fa, fb), (ma, mb)) for (fa, fb), (ma, mb) in zip(*aln.aligned)]
    islands = []
    for blk in blocks:
        if islands:
            (pfa, pfb), (pma, pmb) = islands[-1][-1]
            if blk[0][0] - pfb < min_gap and blk[1][0] - pmb < min_gap:
                islands[-1].append(blk)
                continue
        islands.append([blk])
    f2m = [None] * len(F)
    same = total = 0
    for isl in islands:
        if isl[-1][0][1] - isl[0][0][0] < min_island and len(islands) > 1:
            continue
        for (fa, fb), (ma, mb) in isl:
            for k in range(fb - fa):
                f2m[fa + k] = ma + k
                same += F[fa + k] == M[ma + k]
            total += fb - fa
    return f2m, (same / total if total else 0.0)


def female_specific_segments(f2m, male_len, min_len):
    """Unaligned runs of F of at least `min_len` nt.

    A run at the 5' or 3' end of F counts only if M also has unaligned sequence at that end (a different
    first or last exon); otherwise it is just a longer UTR or contig end and is ignored.
    """
    runs, i, n = [], 0, len(f2m)
    while i < n:
        if f2m[i] is None:
            j = i
            while j < n and f2m[j] is None:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    aligned_m = [m for m in f2m if m is not None]
    if not aligned_m:
        return []
    m_first, m_last = min(aligned_m), max(aligned_m)
    out = []
    for a, b in runs:
        if b - a < min_len:
            continue
        if a == 0 and m_first < min_len:
            continue
        if b == n and male_len - 1 - m_last < min_len:
            continue
        out.append((a, b))
    return out


def cds_from_header(header):
    """1-based CDS start from a 'CDS=start-end' tag in a FASTA header, or None."""
    m = _CDS_TAG.search(header or "")
    return int(m.group(1)) if m else None


def build_pair(pair_id, fid, F, mid, M, male_start=None, female_start=None, min_segment=15, table=1,
               male_header="", female_header=""):
    """Prepare a PairInput from two sequences. Starts are 1-based positions, or None to infer them."""
    notes = []
    f2m, ident = align(F, M)
    if ident < 0.9:
        notes.append(f"aligned identity only {ident:.0%}; check that the two sequences are isoforms of one gene "
                     "and in the same orientation")
    segments = female_specific_segments(f2m, len(M), min_segment)
    if not segments:
        notes.append("no female-specific sequence found")

    if male_start is None:
        male_start = cds_from_header(male_header)
    if male_start is None:
        orf = longest_orf(M, stop_codons(table), table)
        if orf is None:
            raise ValueError(f"{pair_id}: no ATG-initiated ORF with a stop codon in the male sequence {mid}; "
                             "give male_start")
        ms = orf[0]
        notes.append("male start codon = longest ATG ORF of the male sequence")
    else:
        ms = male_start - 1
    if female_start is None:
        female_start = cds_from_header(female_header)
    if female_start is None:
        m2f = {m: f for f, m in enumerate(f2m) if m is not None}
        fs = m2f.get(ms)
        if fs is None or m2f.get(ms + 2) != fs + 2:
            raise ValueError(f"{pair_id}: the male start codon (male position {ms + 1}) is not aligned to the "
                             "female sequence; give female_start")
    else:
        fs = female_start - 1
    return PairInput(pair_id=pair_id, female_id=fid, male_id=mid, female=F, male=M, segments=segments,
                     f2m=f2m, female_start=fs, male_start=ms, mode="transcript", notes=notes)
