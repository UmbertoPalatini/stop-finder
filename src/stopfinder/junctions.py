"""Count, per sample, the assembled sequences that contain each female and male splice junction.

A junction probe (60 nt by default) matches a sequence when it is found on either strand with at most
`max_mismatches` mismatches and no gaps. Candidate positions are found by exact matches of probe blocks
(with k mismatches, at least one of k + 1 blocks matches exactly), then checked base by base.
This records whether each splice form was assembled in each sample; it does not count reads.
"""
from bisect import bisect_right

from .seqio import read_fasta, revcomp


def _hamming_within(a, b, k):
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def contains(probe, seq, k):
    """True if probe matches seq on either strand with at most k mismatches and no gaps."""
    return bool(_matching(probe, seq, [0], k))


def _matching(probe, text, starts, k):
    """Indices of the sequences (concatenated in `text`, sequence i starting at starts[i]) that contain
    probe on either strand with at most k mismatches."""
    n = len(probe)
    size = n // (k + 1)
    hits = set()
    for strand_probe in (probe, revcomp(probe)):
        for off in range(0, (k + 1) * size, size):
            block = strand_probe[off:off + size]
            j = text.find(block)
            while j >= 0:
                s = j - off
                if s >= 0 and "\n" not in text[s:s + n] and _hamming_within(strand_probe, text[s:s + n], k):
                    hits.add(bisect_right(starts, s) - 1)
                j = text.find(block, j + 1)
    return hits


def count(probes, samples, max_mismatches=3):
    """probes: list of {pair_id, junction, sequence}; samples: list of {sample, group, fasta}.

    Returns rows {pair_id, junction, sample, group, n_sequences}.
    """
    rows = []
    for s in samples:
        seqs, _ = read_fasta(s["fasta"])
        starts, pos = [], 0
        for seq in seqs.values():
            starts.append(pos)
            pos += len(seq) + 1
        text = "\n".join(seqs.values())     # one string: a single search covers every sequence
        for p in probes:
            n = len(_matching(p["sequence"], text, starts, max_mismatches))
            rows.append({"pair_id": p["pair_id"], "junction": p["junction"], "sample": s["sample"],
                         "group": s["group"], "n_sequences": n})
    return rows


def summarise(rows):
    """Per pair and junction: number of samples with >= 1 matching sequence, per group."""
    out = {}
    for r in rows:
        key = (r["pair_id"], r["junction"])
        g = out.setdefault(key, {})
        tot, pos = g.get(r["group"], (0, 0))
        g[r["group"]] = (tot + 1, pos + (r["n_sequences"] > 0))
    return [{"pair_id": k[0], "junction": k[1],
             "samples_with_junction": "; ".join(f"{grp}: {pos}/{tot}" for grp, (tot, pos) in sorted(v.items()))}
            for k, v in out.items()]
