#!/usr/bin/env python
"""Extract transcript structures, mRNA, CDS and protein for MOSN (carrie), doublesex and
fruitless from the edited annotations (with the female-specific "FS" transcripts).

Usage: python 01_ref_isoforms.py config.tsv outdir
config.tsv columns: species  gene  gene_id  gtf  genome
Writes: <outdir>/ref_transcripts.tsv, ref_exons.tsv, ref_mrna.fa, ref_cds.fa, ref_protein.fa
"""
import sys, re, subprocess, collections, os
from Bio.Seq import Seq

cfg, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

def faidx(genome, chrom, start, end):
    out = subprocess.run(["samtools", "faidx", genome, f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:]).upper()

def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None

T = open(f"{outdir}/ref_transcripts.tsv", "w")
T.write("species\tgene\tgene_id\ttranscript\tchrom\tstrand\tn_exons\tmrna_len\tcds_len\tprot_len\tcds_complete\n")
E = open(f"{outdir}/ref_exons.tsv", "w")
E.write("species\tgene\ttranscript\texon_rank\tchrom\tstart\tend\tstrand\tlen\n")
FM, FC, FP = (open(f"{outdir}/ref_{x}.fa", "w") for x in ("mrna", "cds", "protein"))

for line in open(cfg):
    if line.startswith("species") or not line.strip():
        continue
    sp, gene, gid, gtf, genome = line.rstrip("\n").split("\t")
    exons, cds = collections.defaultdict(list), collections.defaultdict(list)
    chrom = strand = None
    for l in open(gtf):
        if gid not in l:
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or attr(f[8], "gene_id") != gid:
            continue
        tid = attr(f[8], "transcript_id")
        if tid is None:
            continue
        chrom, strand = f[0], f[6]
        if f[2] == "exon":
            exons[tid].append((int(f[3]), int(f[4])))
        elif f[2] == "CDS":
            cds[tid].append((int(f[3]), int(f[4])))
    lo = min(s for v in exons.values() for s, e in v)
    hi = max(e for v in exons.values() for s, e in v)
    region = faidx(genome, chrom, lo, hi)
    def seg(s, e):
        return region[s - lo:e - lo + 1]
    for tid in sorted(exons):
        ex = sorted(set(exons[tid]), reverse=(strand == "-"))
        def spliced(blocks):
            s = "".join(seg(a, b) for a, b in sorted(blocks))
            return str(Seq(s).reverse_complement()) if strand == "-" else s
        mrna = spliced(ex)
        for i, (a, b) in enumerate(ex, 1):
            E.write(f"{sp}\t{gene}\t{tid}\t{i}\t{chrom}\t{a}\t{b}\t{strand}\t{b - a + 1}\n")
        c = spliced(sorted(set(cds[tid]))) if cds[tid] else ""
        prot = str(Seq(c[: len(c) // 3 * 3]).translate()) if c else ""
        complete = bool(c) and prot.startswith("M") and prot.endswith("*") and "*" not in prot[:-1]
        T.write(f"{sp}\t{gene}\t{gid}\t{tid}\t{chrom}\t{strand}\t{len(ex)}\t{len(mrna)}\t{len(c)}\t{max(len(prot) - 1, 0)}\t{complete}\n")
        name = f"{sp}|{gene}|{tid}"
        FM.write(f">{name}\n{mrna}\n")
        if c:
            FC.write(f">{name}\n{c}\n")
            FP.write(f">{name}\n{prot}\n")
for fh in (T, E, FM, FC, FP):
    fh.close()
