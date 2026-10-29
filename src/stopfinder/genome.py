"""Genome mode: female and male mRNAs are built from a genome FASTA and a GTF/GFF3 annotation.

The female-specific sequence is the set of genomic positions in the female mRNA that are absent from the
male mRNA.
"""
from .core import PairInput
from .seqio import revcomp


def merge(blocks):
    """Merge overlapping or touching exons (keeps the longest span)."""
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def parse_exons(text, chrom):
    """'start-end;start-end' or 'chrom:start-end;...' -> [(start, end)]. '.' or '' -> []."""
    out = []
    for part in (text or "").replace(",", ";").split(";"):
        part = part.strip()
        if not part or part == ".":
            continue
        if ":" in part:
            c, part = part.rsplit(":", 1)
            if c != chrom:
                raise ValueError(f"extra exon {c}:{part} is not on the transcript's chromosome {chrom}")
        s, e = (int(x.replace(",", "")) for x in part.split("-"))
        out.append((min(s, e), max(s, e)))
    return out


def spliced(fasta, chrom, strand, blocks):
    """mRNA sequence and the genomic coordinate of each base, both in transcript orientation."""
    seq = "".join(fasta.fetch(chrom, s, e) for s, e in blocks)
    pos = [p for s, e in blocks for p in range(s, e + 1)]
    if len(seq) != len(pos):
        raise ValueError(f"exons extend beyond the end of {chrom}")
    if strand == "-":
        return revcomp(seq), pos[::-1]
    return seq, pos


def build_pair(pair_id, fasta, ftx, mtx, extra_exons="", start=None):
    """ftx, mtx: annotation.Transcript. start: genomic coordinate of the shared start codon, or None
    to use the annotated CDS start of the male transcript."""
    notes = []
    if (ftx.chrom, ftx.strand) != (mtx.chrom, mtx.strand):
        raise ValueError(f"{pair_id}: {ftx.id} and {mtx.id} are on different chromosomes or strands")
    chrom, strand = mtx.chrom, mtx.strand
    extra = parse_exons(extra_exons, chrom)
    fb = merge(ftx.exons + extra)
    mb = merge(mtx.exons)
    if extra:
        notes.append(f"female exons added: {extra_exons}")
    F, fpos = spliced(fasta, chrom, strand, fb)
    M, mpos = spliced(fasta, chrom, strand, mb)
    if start is None:
        start = mtx.cds_start()
        if start is None:
            raise ValueError(f"{pair_id}: male transcript {mtx.id} has no CDS in the annotation; give start")
    m_index = {p: i for i, p in enumerate(mpos)}
    try:
        fs, ms = fpos.index(start), m_index[start]
    except (ValueError, KeyError):
        raise ValueError(f"{pair_id}: start codon {chrom}:{start} is not in both transcripts") from None
    f2m = [m_index.get(p) for p in fpos]
    segments, i = [], 0
    while i < len(F):
        if f2m[i] is None:
            j = i + 1
            while j < len(F) and f2m[j] is None and abs(fpos[j] - fpos[j - 1]) == 1:
                j += 1
            segments.append((i, j))
            i = j
        else:
            i += 1

    def coord(i):
        return f"{chrom}:{fpos[i]}({strand})"

    def seg_coord(a, b):
        lo, hi = min(fpos[a], fpos[b - 1]), max(fpos[a], fpos[b - 1])
        return f"{chrom}:{lo}-{hi}"

    return PairInput(pair_id=pair_id, female_id=ftx.id, male_id=mtx.id, female=F, male=M, segments=segments,
                     f2m=f2m, female_start=fs, male_start=ms, mode="genome", female_coord=coord,
                     segment_coord=seg_coord, notes=notes)
