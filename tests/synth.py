"""Synthetic genes with known answers, used by the tests and by examples/make_example.py."""
import random

from stopfinder.seqio import revcomp

STOPS = {"TAA", "TAG", "TGA"}
_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def coding(n_codons, rng):
    """Random codons without stops."""
    return "".join(rng.choice(_SENSE) for _ in range(n_codons))


def stopless(n, rng):
    """Random sequence of length n with no stop codon (and no ATG) in any frame."""
    s = rng.choice("ACG")
    while len(s) < n:
        c = rng.choice("ACGT")
        if s[-2:] + c not in STOPS | {"ATG"}:
            s += c
    return s


def noise(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def gene_model(seed=1):
    """Exon sequences of a gene with a female-specific exon.

    Male:    E1 E2 E3 E4         (CDS starts in E1, stop in E4)
    Female:  E1 E2 Ef E3 E4      Ef = female exon; variants below
    All exon lengths are multiples of 3 except where stated, and every exon is in frame 0 of the ORF.
    """
    rng = random.Random(seed)
    utr5 = stopless(40, rng)
    e1 = utr5 + "ATG" + coding(40, rng)              # 163 nt
    e2 = coding(50, rng)                             # 150
    e3 = coding(60, rng)                             # 180
    e4 = coding(30, rng) + "TAA" + noise(80, rng)    # male stop + 3' UTR
    # female exon A: 20 codons, then TGA, TAG (double stop, as in dsx), then 10 codons
    efa = coding(20, rng) + "TGATAG" + coding(10, rng)   # 96 nt
    # female exon B: 220 nt (not a multiple of 3), no stop in any frame -> frameshift
    efb = stopless(220, rng)
    return dict(utr5=len(utr5), e1=e1, e2=e2, e3=e3, e4=e4, efa=efa, efb=efb, rng=rng)


def build_genome(m, seed=2):
    """Place the gene on the + strand of chr1 and its reverse complement on the - strand of chr2.

    Returns genome FASTA text, GTF text, and the exon coordinates.
    """
    rng = random.Random(seed)
    order = ["e1", "e2", "efa", "efb", "e3", "e4"]
    pos, seq, coords = 1, [], {}
    intron = lambda: "GT" + noise(196, rng) + "AG"
    lead = noise(500, rng)
    seq.append(lead)
    pos += len(lead)
    for k in order:
        coords[k] = (pos, pos + len(m[k]) - 1)
        seq.append(m[k])
        pos += len(m[k])
        it = intron()
        seq.append(it)
        pos += len(it)
    chr1 = "".join(seq) + noise(500, rng)
    chr2 = revcomp(chr1)
    L = len(chr1)
    gtf = []

    def tx(chrom, strand, tid, exons, cds_start, cds_end):
        for k in exons:
            s, e = coords[k]
            if strand == "-":
                s, e = L - e + 1, L - s + 1
            gtf.append(f'{chrom}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "G"; transcript_id "{tid}";')
        if cds_start:
            s, e = cds_start, cds_end
            if strand == "-":
                s, e = L - e + 1, L - s + 1
            gtf.append(f'{chrom}\ttest\tCDS\t{s}\t{e}\t.\t{strand}\t0\tgene_id "G"; transcript_id "{tid}";')

    atg = coords["e1"][0] + m["utr5"]
    male_stop_end = coords["e4"][0] + 30 * 3 + 2
    for chrom, strand, sfx in (("chr1", "+", ""), ("chr2", "-", ".minus")):
        # CDS given as one span: only its first base (the start codon) is used
        tx(chrom, strand, "M" + sfx, ["e1", "e2", "e3", "e4"], atg, male_stop_end)
        tx(chrom, strand, "FA" + sfx, ["e1", "e2", "efa", "e3", "e4"], None, None)
        tx(chrom, strand, "FB" + sfx, ["e1", "e2", "efb", "e3", "e4"], None, None)
    fasta = f">chr1\n" + "\n".join(chr1[i:i + 60] for i in range(0, L, 60)) + \
            f"\n>chr2\n" + "\n".join(chr2[i:i + 60] for i in range(0, L, 60)) + "\n"
    return fasta, "\n".join(gtf) + "\n", coords, L


def mrnas(m):
    male = m["e1"] + m["e2"] + m["e3"] + m["e4"]
    fa = m["e1"] + m["e2"] + m["efa"] + m["e3"] + m["e4"]
    fb = m["e1"] + m["e2"] + m["efb"] + m["e3"] + m["e4"]
    return male, fa, fb


def random_pair(seed, snp_every=50):
    """A random female/male pair of one of three structures, with SNPs in the shared 5' part of F.

    kind 0: cassette female exon containing a stop; kind 1: alternative last exon (female ends in its own
    exon); kind 2: mutually exclusive exons. Returns kind, F, M, and the female-specific interval (a, b).
    """
    rng = random.Random(seed)
    shared = stopless(50, rng) + "ATG" + coding(rng.randint(30, 300), rng)
    down = coding(rng.randint(30, 200), rng) + "TAA" + noise(rng.randint(20, 300), rng)
    kind = seed % 3
    if kind == 0:
        fe = coding(rng.randint(5, 60), rng) + rng.choice(sorted(STOPS)) + noise(rng.randint(10, 150), rng)
        F, M = shared + fe + down, shared + down
    elif kind == 1:
        fe = coding(rng.randint(5, 60), rng) + "TGA" + noise(rng.randint(50, 300), rng)
        F, M = shared + fe, shared + coding(80, rng) + "TAA" + noise(100, rng)
    else:
        fe = coding(rng.randint(5, 40), rng) + "TAG" + coding(rng.randint(5, 40), rng)
        F, M = shared + fe + down, shared + coding(rng.randint(20, 60), rng) + down
    f = list(F)
    n = 0
    while n < max(3, len(shared) // snp_every):
        i = rng.randrange(60, len(shared) - 5)
        old = f[i]
        f[i] = rng.choice([c for c in "ACGT" if c != old])
        c0 = 53 + (i - 53) // 3 * 3            # codon of the ORF (ATG at 50) containing position i
        if "".join(f[c0:c0 + 3]) in STOPS:      # do not create a premature stop
            f[i] = old
            continue
        n += 1
    return kind, "".join(f), M, (len(shared), len(shared) + len(fe))
