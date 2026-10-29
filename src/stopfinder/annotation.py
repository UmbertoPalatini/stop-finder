"""Read exon and CDS structure of selected transcripts from a GTF or GFF3 file."""
import gzip
import re
from dataclasses import dataclass, field

_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')
# Prefixes that GFF3 files from NCBI and Ensembl put in front of transcript IDs.
_PREFIXES = ("transcript:", "rna-", "mRNA:", "rna:")


@dataclass
class Transcript:
    id: str
    chrom: str = None
    strand: str = None
    exons: list = field(default_factory=list)   # (start, end), 1-based inclusive
    cds: list = field(default_factory=list)

    def cds_start(self):
        """Genomic coordinate of the first base of the start codon, or None."""
        if not self.cds:
            return None
        return min(s for s, _ in self.cds) if self.strand == "+" else max(e for _, e in self.cds)


def _aliases(tid):
    out = {tid}
    for p in _PREFIXES:
        if tid.startswith(p):
            out.add(tid[len(p):])
    return out


def _transcript_ids(attrs):
    """Transcript IDs of a feature line (GTF: transcript_id; GFF3: Parent, possibly several)."""
    if '"' in attrs or "transcript_id " in attrs:
        m = dict(_GTF_ATTR.findall(attrs))
        return [m["transcript_id"]] if "transcript_id" in m else []
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith("Parent="):
            return part[7:].split(",")
        if part.startswith("transcript_id="):
            return [part[14:]]
    return []


def load_transcripts(path, wanted):
    """Return {wanted_id: Transcript} for the requested transcript IDs.

    IDs are matched exactly, or after removing common GFF3 prefixes such as 'rna-' or 'transcript:'.
    """
    wanted = set(wanted)
    found = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            for raw in _transcript_ids(f[8]):
                hit = wanted & _aliases(raw)
                for tid in hit:
                    t = found.setdefault(tid, Transcript(tid))
                    if t.chrom is not None and (t.chrom, t.strand) != (f[0], f[6]):
                        raise ValueError(f"{path}: transcript {tid} has features on more than one chromosome/strand")
                    t.chrom, t.strand = f[0], f[6]
                    (t.exons if f[2] == "exon" else t.cds).append((int(f[3]), int(f[4])))
    for t in found.values():
        if not t.exons:          # some files list only CDS features
            t.exons = list(t.cds)
        t.exons = sorted(set(t.exons))
        t.cds = sorted(set(t.cds))
    return found
