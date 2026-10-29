#!/usr/bin/env python
"""Female vs male isoforms and stop codons in the female-specific exonic sequence,
for species with a reference genome (Aedes, Culex, Anopheles).

Usage: python 02_ref_stops.py pairs.tsv ref_genes.tsv ref_exons.tsv outdir
pairs.tsv columns: species gene female_tx male_tx add_exons(chrom-free 'start-end;...' or '.')
  female chain = exons of female_tx + add_exons, overlapping exons merged (longest span kept).
  ORF start   = first base of the annotated CDS of male_tx (shared N-terminus).
Writes <outdir>/ref_pairs.tsv, ref_female_specific_stops.tsv, ref_FM_mrna.fa, ref_FM_protein.fa
"""
import sys, re, subprocess, collections, os
from Bio.Seq import Seq

pairs_f, genes_f, exons_f, outdir = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)
STOPS = {"TAA", "TAG", "TGA"}

genes = {}
for l in open(genes_f):
    f = l.rstrip("\n").split("\t")
    if f[0] != "species":
        genes[(f[0], f[1])] = f
ex = collections.defaultdict(list)
meta = {}
for l in open(exons_f):
    f = l.rstrip("\n").split("\t")
    if f[0] == "species":
        continue
    ex[(f[0], f[1], f[2])].append((int(f[5]), int(f[6])))
    meta[(f[0], f[1])] = (f[4], f[7])

def faidx(genome, chrom, start, end):
    out = subprocess.run(["samtools", "faidx", genome, f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:]).upper()

def cds_start(gtf, gid, tid, strand):
    pos = []
    for l in open(gtf):
        if tid not in l:
            continue
        f = l.split("\t")
        if f[2] == "CDS" and f'transcript_id "{tid}"' in f[8]:
            pos += [int(f[3]), int(f[4])]
    return min(pos) if strand == "+" else max(pos)

def merge(blocks):
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out

def chain_positions(blocks, strand):
    """genomic coordinate of each mRNA base, in transcript orientation"""
    pos = [p for s, e in sorted(blocks) for p in range(s, e + 1)]
    return pos[::-1] if strand == "-" else pos

P = open(f"{outdir}/ref_pairs.tsv", "w")
P.write("species\tgene\tfemale_chain\tmale_tx\tchrom\tstrand\tstart_codon_pos\tfemale_prot_len\tmale_prot_len\t"
        "female_specific_segments\tfemale_stop_codon\tfemale_stop_genomic\tstop_in_female_specific\tstop_context_-6_+6\t"
        "nt4_after_stop\tcodons_to_next_inframe_stop\tnext_stop_codon\tfemale_specific_inframe_stops\n")
S = open(f"{outdir}/ref_female_specific_stops.tsv", "w")
S.write("species\tgene\tsegment\tframe_rel_ORF\tstop_index\tcodon\tgenomic_pos\tcontext_-6_+6\n")
FM = open(f"{outdir}/ref_FM_mrna.fa", "w")
FP = open(f"{outdir}/ref_FM_protein.fa", "w")

for l in open(pairs_f):
    f = l.rstrip("\n").split("\t")
    if f[0] == "species" or not l.strip():
        continue
    sp, gene, ftx, mtx, add = f[:5]
    chrom, strand = meta[(sp, gene)]
    gtf, genome = genes[(sp, gene)][3], genes[(sp, gene)][4]
    gid = genes[(sp, gene)][2]
    fb = list(ex[(sp, gene, ftx)])
    if add != ".":
        fb += [tuple(map(int, x.split("-"))) for x in add.split(";")]
    fb = merge(fb)
    mb = merge(ex[(sp, gene, mtx)])
    lo = min(s for s, e in fb + mb); hi = max(e for s, e in fb + mb)
    region = faidx(genome, chrom, lo, hi)
    def mrna(blocks):
        pos = chain_positions(blocks, strand)
        seq = "".join(region[p - lo] for p in pos)
        if strand == "-":
            seq = str(Seq(seq).complement())
        return seq, pos
    fseq, fpos = mrna(fb)
    mseq, mpos = mrna(mb)
    st = cds_start(gtf, gid, mtx, strand)
    fi, mi = fpos.index(st), mpos.index(st)
    def orf(seq, i):
        cod = [seq[j:j + 3] for j in range(i, len(seq) - 2, 3)]
        for k, c in enumerate(cod):
            if c in STOPS:
                return k, c
        return None, None
    fk, fstop = orf(fseq, fi)
    mk, mstop = orf(mseq, mi)
    fprot = str(Seq(fseq[fi:fi + 3 * (fk + 1)]).translate()) if fk is not None else ""
    mprot = str(Seq(mseq[mi:mi + 3 * (mk + 1)]).translate()) if mk is not None else ""
    mset = set(mpos)
    fspec = [j for j, p in enumerate(fpos) if p not in mset]
    # contiguous female-specific segments in mRNA coordinates
    segs = []
    for j in fspec:
        if segs and j == segs[-1][1] + 1 and abs(fpos[j] - fpos[j - 1]) == 1:
            segs[-1][1] = j
        else:
            segs.append([j, j])
    seg_txt = ";".join(f"{chrom}:{min(fpos[a], fpos[b])}-{max(fpos[a], fpos[b])}({b - a + 1}bp)" for a, b in segs)
    fstop_idx = fi + 3 * fk if fk is not None else None
    in_spec = fstop_idx is not None and all((fstop_idx + d) in set(fspec) for d in range(3))
    ctx = fseq[fstop_idx - 6: fstop_idx + 9] if fstop_idx else ""
    nt4 = fseq[fstop_idx + 3] if fstop_idx and fstop_idx + 3 < len(fseq) else ""
    # next in-frame stop after the first one (read-through extension)
    nxt = nxtc = ""
    if fstop_idx is not None:
        for j in range(fstop_idx + 3, len(fseq) - 2, 3):
            if fseq[j:j + 3] in STOPS:
                nxt, nxtc = (j - fstop_idx) // 3, fseq[j:j + 3]
                break
    # all in-frame (ORF frame) stops inside each female-specific segment, plus stops in the other 2 frames
    inframe = []
    for si, (a, b) in enumerate(segs, 1):
        for j in range(a, b - 1):
            c = fseq[j:j + 3]
            if c in STOPS and j + 2 <= b:
                frame = (j - fi) % 3
                if frame == 0:
                    inframe.append(c)
                S.write(f"{sp}\t{gene}\t{si}\t{frame}\t{(j - fi) // 3 if frame == 0 else ''}\t{c}\t{fpos[j]}\t{fseq[max(j - 6, 0):j + 9]}\n")
    cnt = collections.Counter(inframe)
    P.write(f"{sp}\t{gene}\t{ftx}{'+' + add if add != '.' else ''}\t{mtx}\t{chrom}\t{strand}\t{st}\t{len(fprot) - 1}\t{len(mprot) - 1}\t"
            f"{seg_txt}\t{fstop}\t{fpos[fstop_idx] if fstop_idx is not None else ''}\t{in_spec}\t{ctx}\t{nt4}\t{nxt}\t{nxtc}\t"
            f"{','.join(f'{k}:{v}' for k, v in sorted(cnt.items()))}\n")
    FM.write(f">{sp}|{gene}|female|{ftx}\n{fseq}\n>{sp}|{gene}|male|{mtx}\n{mseq}\n")
    FP.write(f">{sp}|{gene}|female\n{fprot}\n>{sp}|{gene}|male\n{mprot}\n")
for h in (P, S, FM, FP):
    h.close()
