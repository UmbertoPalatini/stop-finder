"""FASTA reading and random access to indexed genomes (samtools-compatible .fai)."""
import gzip
import os

_COMP = str.maketrans("ACGTUNacgtunRYKMSWBDHVrykmswbdhv", "TGCAANtgcaanYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """Return {id: sequence} (id = first word of the header) and {id: full header}."""
    seqs, headers, name, buf = {}, {}, None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                header = line[1:].strip()
                name = header.split()[0] if header else ""
                if name in seqs:
                    raise ValueError(f"{path}: duplicate sequence id '{name}'")
                headers[name], buf = header, []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs, headers


class IndexedFasta:
    """Random access to an uncompressed FASTA through a samtools-style .fai index.

    The index is built (and saved next to the FASTA, if writable) when it does not exist.
    """

    def __init__(self, path):
        if path.endswith(".gz"):
            raise ValueError(f"{path}: genome FASTA must be uncompressed (or bgzip + samtools faidx; "
                             "plain gzip cannot be indexed)")
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai) or os.path.getmtime(fai) < os.path.getmtime(path):
            self.index = self._build()
            try:
                with open(fai, "w") as out:
                    for name, (length, offset, blen, llen) in self.index.items():
                        out.write(f"{name}\t{length}\t{offset}\t{blen}\t{llen}\n")
            except OSError:
                pass
        else:
            self.index = {}
            with open(fai) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    self.index[f[0]] = tuple(int(x) for x in f[1:5])
        self.fh = open(path, "rb")

    def _build(self):
        index, name, length, offset, blen, llen = {}, None, 0, 0, None, None
        pos, short_line = 0, False
        with open(self.path, "rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        index[name] = (length, offset, blen or 0, llen or 0)
                    name = raw[1:].split()[0].decode()
                    length, blen, llen, offset, short_line = 0, None, None, pos + len(raw), False
                else:
                    seq = raw.rstrip(b"\r\n")
                    if blen is None:
                        blen, llen = len(seq), len(raw)
                    elif (short_line and seq) or len(seq) > blen:
                        raise ValueError(f"{self.path}: sequence '{name}' has lines of different lengths; "
                                         "reformat it (e.g. seqkit seq -w 60) before use")
                    if len(seq) < blen:
                        short_line = True       # only the last line of a sequence may be shorter
                    length += len(seq)
                pos += len(raw)
        if name is not None:
            index[name] = (length, offset, blen or 0, llen or 0)
        return index

    def fetch(self, chrom, start, end):
        """Sequence of chrom:start-end (1-based, inclusive), upper case."""
        if chrom not in self.index:
            raise KeyError(f"sequence '{chrom}' not found in {self.path}")
        length, offset, blen, llen = self.index[chrom]
        start, end = max(start, 1), min(end, length)
        if end < start:
            return ""
        s0, e0 = start - 1, end
        b0 = offset + (s0 // blen) * llen + s0 % blen
        b1 = offset + ((e0 - 1) // blen) * llen + (e0 - 1) % blen + 1
        self.fh.seek(b0)
        data = self.fh.read(b1 - b0).decode()
        return data.replace("\n", "").replace("\r", "").upper()

    def close(self):
        self.fh.close()
