#!/usr/bin/env python
"""Find MOSN, dsx and fru contigs in the per-sample Trinity assemblies of Toxorhynchites and Wyeomyia
and classify each as female-type or male-type at the sex-specific splice junction.

Logic: the female and male reference proteins share an N-terminus of `shared` aa (from 02). For each
candidate contig, align the reference male proteins (tblastn, with sequences) and locate the contig codon
aligned to the last shared residue (J). If the male alignment continues >= MIN_CONT aa past J in the same
HSP, the contig is male-type. Otherwise read the contig in frame from J to the first stop: if this happens
before the male alignment resumes, the contig is female-type and the peptide after J is its female-specific tail.

Usage: python 04_classify_contigs.py samples.tsv blastdir dbdir ref_FM_protein.fa outdir
"""
import sys, os, re, subprocess, collections
from Bio.Seq import Seq

samples_f, blastdir, dbdir, protfa, outdir = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)
MIN_SCORE, MIN_CONT, STOPS = 80, 15, {"TAA", "TAG", "TGA"}

def readfa(path):
    d, h = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith(">"):
            h = l[1:].split()[0]; d[h] = []
        elif h:
            d[h].append(l)
    return {k: "".join(v) for k, v in d.items()}

prot = readfa(protfa)
shared = {}
for h, s in prot.items():
    sp, g, sex = h.split("|")
    if sex == "female":
        m = prot[f"{sp}|{g}|male"]
        i = 0
        while s[i] == m[i]:
            i += 1
        shared[(sp, g)] = i          # number of shared N-terminal residues
genes = sorted({g for _, g in shared})

samples = [l.rstrip("\n").split("\t") for l in open(samples_f) if not l.startswith("species")]
OUT = open(f"{outdir}/contig_classes.tsv", "w")
OUT.write("species\tsample\tsex\tgene\tcontig\tcontig_len\tscore\tbest_ref\tstrand\tref_cov_male\tclass\tJ_contig\t"
          "aa_after_J_to_stop\tstop_codon\tstop_context_-6_+6\tnext_inframe_stop_codons\ttail_peptide\torf_len_aa\torf_complete\torf_protein\n")
FA = open(f"{outdir}/candidate_contigs_oriented.fa", "w")

def tblastn_subject(query_fa, subj_fa):
    out = subprocess.run(["tblastn", "-query", query_fa, "-subject", subj_fa, "-evalue", "1e-3", "-outfmt",
                          "6 qseqid sseqid qstart qend sstart send bitscore sframe qseq sseq"],
                         capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.strip().split("\n") if l]

tmpq, tmps = f"{outdir}/.q.fa", f"{outdir}/.s.fa"
for species, sample, sex in samples:
    hits = collections.defaultdict(lambda: collections.defaultdict(float))
    for l in open(f"{blastdir}/{sample}.tblastn.tsv"):
        f = l.rstrip("\n").split("\t")
        q, s, bits = f[0], f[1], float(f[9])
        hits[s][q] += bits
    best = {}
    for s, qd in hits.items():
        q, sc = max(qd.items(), key=lambda x: x[1])
        if sc >= MIN_SCORE:
            best[s] = (q.split("|")[1], q, sc)
    if not best:
        continue
    ids = f"{outdir}/.ids"
    open(ids, "w").write("\n".join(best) + "\n")
    seqs = subprocess.run(["blastdbcmd", "-db", f"{dbdir}/{sample}", "-entry_batch", ids],
                          capture_output=True, text=True, check=True).stdout
    cseq = {}
    h = None
    for l in seqs.split("\n"):
        if l.startswith(">"):
            h = l[1:].split()[0]; cseq[h] = ""
        elif h:
            cseq[h] += l.strip()
    for c, (g, bq, sc) in sorted(best.items(), key=lambda x: (x[1][0], -x[1][2])):
        seq = cseq.get(c) or cseq.get(c.split("|")[-1])
        if seq is None:
            continue
        # align the three reference male proteins of this gene to the contig, keep the best
        with open(tmpq, "w") as fq:
            for (sp, gg) in shared:
                if gg == g:
                    fq.write(f">{sp}|{g}|male\n{prot[f'{sp}|{g}|male'].rstrip('*')}\n")
        open(tmps, "w").write(f">c\n{seq}\n")
        hs = tblastn_subject(tmpq, tmps)
        if not hs:
            continue
        tot = collections.defaultdict(float)
        for x in hs:
            tot[x[0]] += float(x[6])
        ref = max(tot, key=tot.get)
        hs = [x for x in hs if x[0] == ref]
        strand = "+" if int(max(hs, key=lambda x: float(x[6]))[7]) > 0 else "-"
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
            open(tmps, "w").write(f">c\n{seq}\n")
            open(tmpq, "w").write(f">{ref}\n{prot[ref].rstrip('*')}\n")
            hs = [x for x in tblastn_subject(tmpq, tmps) if int(x[7]) > 0]
        else:
            hs = [x for x in hs if int(x[7]) > 0]
        FA.write(f">{species}|{sample}|{sex}|{g}|{c}\n{seq}\n")
        sh = shared[tuple(ref.split("|")[:2])]
        qcov = sorted({p for x in hs for p in range(int(x[2]), int(x[3]) + 1)})
        cov = f"{qcov[0]}-{qcov[-1]}({len(qcov)}/{len(prot[ref]) - 1})" if qcov else ""
        # map query residue `sh` (last shared) to contig codon, and test continuation
        J = cont = None
        for x in sorted(hs, key=lambda x: -float(x[6])):
            qs, qe, ss = int(x[2]), int(x[3]), int(x[4])
            if not (qs <= sh <= qe + 3 and qe >= sh - 40):
                continue
            qi, si = qs, ss
            last = None
            for a, b in zip(x[8], x[9]):
                if a != "-":
                    if qi <= sh:
                        last = si if b != "-" else last
                if b != "-":
                    si += 3
                if a != "-":
                    qi += 1
            if last is None:
                continue
            J = last + 2                      # last nt of codon aligned to residue `sh` (or closest before)
            cont = qe - sh
            break
        cls, tail, stopc, ctx, nxt, after = "no_junction", "", "", "", "", ""
        if J is not None:
            if cont >= MIN_CONT:
                cls = "male"
            else:
                codons = []      # J (1-based last nt of junction codon) == 0-based index of the next nt
                j = J
                while j + 3 <= len(seq):
                    cdn = seq[j:j + 3]
                    codons.append(cdn)
                    if cdn in STOPS:
                        break
                    j += 3
                if codons and codons[-1] in STOPS:
                    cls = "female"
                    tail = str(Seq("".join(codons)).translate())
                    stopc, after = codons[-1], len(codons) - 1
                    ctx = seq[max(j - 6, 0):j + 9]
                    nx = [seq[t:t + 3] for t in range(j + 3, min(len(seq) - 2, j + 3 + 3 * 300), 3) if seq[t:t + 3] in STOPS]
                    nxt = ",".join(nx[:3])
                else:
                    cls = "female_or_truncated(no_stop)"
                    tail = str(Seq("".join(codons)[: len(codons) * 3 // 3 * 3]).translate())
        # ORF containing the aligned region (frame of best HSP)
        bh = max(hs, key=lambda x: float(x[6]))
        fr = (int(bh[4]) - 1) % 3
        aa = str(Seq(seq[fr: fr + (len(seq) - fr) // 3 * 3]).translate())
        mid = (int(bh[4]) - 1 - fr) // 3
        a0 = aa.rfind("*", 0, mid) + 1
        b0 = aa.find("*", mid)
        orf = aa[a0: b0 + 1 if b0 >= 0 else len(aa)]
        m = orf.find("M")
        orf_m = orf[m:] if m >= 0 else orf
        complete = orf_m.startswith("M") and orf_m.endswith("*") and a0 > 0
        OUT.write(f"{species}\t{sample}\t{sex}\t{g}\t{c}\t{len(seq)}\t{sc:.0f}\t{ref}\t{strand}\t{cov}\t{cls}\t{J or ''}\t"
                  f"{after}\t{stopc}\t{ctx}\t{nxt}\t{tail}\t{len(orf_m.rstrip('*'))}\t{complete}\t{orf_m}\n")
OUT.close(); FA.close()
for t in (tmpq, tmps, f"{outdir}/.ids"):
    if os.path.exists(t):
        os.remove(t)
