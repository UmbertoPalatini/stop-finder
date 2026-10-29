#!/usr/bin/env python
"""Female/male transcripts and female-specific stop codons for Toxorhynchites and Wyeomyia.

For each representative pair (female contig F, male contig M; chosen from 04):
  1. blastn F vs M. The female-specific segment is the part of F between two HSPs that are adjacent on M
     (an insertion in F at a male splice junction).
  2. The male ORF is the longest ATG-initiated ORF of M. Its start is mapped onto F, and F is translated
     from there to the first stop (female protein).
  3. Stop codon, sequence context, next in-frame stop, and all in-frame stops in the female segment.
  4. Sex support: the female junction (60 nt: 30 shared + 30 female segment) and the male junction
     (30 nt either side of the male splice site) are searched (blastn-short, >=95% id over >=56 nt)
     in every candidate contig of that gene from all 12 assemblies (from 04).

Usage: python 05_toxwye_pairs.py pairs.tsv candidate_contigs_oriented.fa samples.tsv outdir
pairs.tsv: species gene isoform female_contig_header male_contig_header   (headers as in the 04 fasta)
"""
import sys, os, subprocess, collections, tempfile
from Bio.Seq import Seq

pairs_f, cand_f, samples_f, outdir = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)
STOPS = {"TAA", "TAG", "TGA"}

def readfa(p):
    d, h = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            h = l[1:]; d[h] = []
        elif h:
            d[h].append(l)
    return {k: "".join(v) for k, v in d.items()}

cand = readfa(cand_f)
sex_of = {l.split("\t")[1]: l.rstrip("\n").split("\t")[2] for l in open(samples_f) if not l.startswith("species")}
tmp = tempfile.mkdtemp(dir=outdir)

def blastn(q, s, short=False):
    open(f"{tmp}/q.fa", "w").write(f">q\n{q}\n")
    open(f"{tmp}/s.fa", "w").write(s if s.startswith(">") else f">s\n{s}\n")
    cmd = ["blastn", "-query", f"{tmp}/q.fa", "-subject", f"{tmp}/s.fa", "-outfmt",
           "6 qstart qend sstart send pident length sseqid", "-evalue", "1e-5" if short else "1e-10", "-strand", "plus"]
    if short:
        cmd += ["-task", "blastn-short", "-word_size", "11"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.strip().split("\n") if l]

def longest_orf(s):
    best = (0, None, None)
    for fr in range(3):
        aa = str(Seq(s[fr:fr + (len(s) - fr) // 3 * 3]).translate())
        i = 0
        for seg in aa.split("*"):
            m = seg.find("M")
            if m >= 0 and len(seg) - m > best[0] and i + len(seg) < len(aa):   # require a stop codon
                best = (len(seg) - m, fr + 3 * (i + m), fr + 3 * (i + len(seg)) + 3)
            i += len(seg) + 1
    return best   # aa length, 0-based start, 0-based end (exclusive, includes stop)

P = open(f"{outdir}/toxwye_pairs.tsv", "w")
P.write("species\tgene\tisoform\tfemale_contig\tmale_contig\tfemale_prot_len\tmale_prot_len\tfemale_specific_segments_in_F\t"
        "female_stop_codon\tfemale_stop_pos_in_F\tstop_in_female_specific\tstop_context_-6_+6\tnt4_after_stop\t"
        "codons_to_next_inframe_stop\tnext_stop_codon\tfemale_specific_inframe_stops\tfemale_tail_peptide\t"
        "female_junction_samples(F/M)\tmale_junction_samples(F/M)\n")
S = open(f"{outdir}/toxwye_female_specific_stops.tsv", "w")
S.write("species\tgene\tisoform\tsegment\tframe_rel_ORF\tstop_index\tcodon\tpos_in_F\tcontext_-6_+6\n")
J = open(f"{outdir}/toxwye_junction_support.tsv", "w")
J.write("species\tgene\tisoform\tjunction\tsample\tsex\tn_contigs\n")
TM = open(f"{outdir}/toxwye_transcripts.fa", "w")
TP = open(f"{outdir}/toxwye_proteins.fa", "w")

for l in open(pairs_f):
    f = l.rstrip("\n").split("\t")
    if f[0] == "species" or not l.strip():
        continue
    sp, gene, iso, fh, mh = f[:5]
    F, M = cand[fh], cand[mh]
    hs = sorted(([int(x[0]), int(x[1]), int(x[2]), int(x[3]), float(x[4])] for x in blastn(F, M)), key=lambda x: x[0])
    hs = [h for h in hs if h[2] < h[3]]
    # collinear chain on both sequences (greedy by query start)
    chain = []
    for h in hs:
        if not chain or (h[0] > chain[-1][0] and h[2] > chain[-1][2]):
            chain.append(h)
    segs = []   # female-specific segments: gaps in F where M is (nearly) contiguous
    for a, b in zip(chain, chain[1:]):
        gapF, gapM = b[0] - a[1] - 1, b[2] - a[3] - 1
        if gapF > 15 and abs(gapM) <= 15:
            segs.append((a[1] + 1, b[0] - 1, a[3]))     # 1-based F coords, M junction (last M base before)
    if not segs and chain and chain[-1][1] < len(F) - 30:   # female sequence diverges to the 3' end
        segs.append((chain[-1][1] + 1, len(F), chain[-1][3]))
    # male ORF and its start on F
    mL, ms, me = longest_orf(M)
    mprot0 = str(Seq(M[ms:me]).translate())
    fs = None                                   # male start on F: find the male N-terminal peptide in F
    for off in range(0, 40, 4):
        pep = mprot0[off:off + 12]
        for fr in range(3):
            aa = str(Seq(F[fr:fr + (len(F) - fr) // 3 * 3]).translate())
            i = aa.find(pep)
            if i >= 0:
                fs = fr + 3 * (i - off)
                break
        if fs is not None:
            break
    cod, k = [], fs
    while k is not None and k + 3 <= len(F):
        cod.append(F[k:k + 3])
        if F[k:k + 3] in STOPS:
            break
        k += 3
    stop_i = k if cod and cod[-1] in STOPS else None
    fprot = str(Seq("".join(cod)).translate()) if stop_i is not None else ""
    spec_pos = {p for a, b, _ in segs for p in range(a - 1, b)}
    mprot = mprot0
    in_spec = stop_i is not None and all(p in spec_pos for p in range(stop_i, stop_i + 3))
    ctx = F[max(stop_i - 6, 0):stop_i + 9] if stop_i else ""
    nt4 = F[stop_i + 3] if stop_i and stop_i + 3 < len(F) else ""
    nxt = nxtc = ""
    if stop_i is not None:
        for j in range(stop_i + 3, len(F) - 2, 3):
            if F[j:j + 3] in STOPS:
                nxt, nxtc = (j - stop_i) // 3, F[j:j + 3]
                break
    inframe = []
    for si, (a, b, _) in enumerate(segs, 1):
        for j in range(a - 1, b - 2):
            c = F[j:j + 3]
            if c in STOPS:
                fr = (j - fs) % 3 if fs is not None else ""
                if fr == 0:
                    inframe.append(c)
                S.write(f"{sp}\t{gene}\t{iso}\t{si}\t{fr}\t{(j - fs) // 3 if fr == 0 else ''}\t{c}\t{j + 1}\t{F[max(j - 6, 0):j + 9]}\n")
    # the female-specific segment that carries the stop (others are listed as other differences)
    key = None
    if stop_i is not None:
        before = [s_ for s_ in segs if s_[0] - 1 <= stop_i]
        key = before[-1] if before else None
    tail = ""
    if key and fs is not None:
        first = key[0] - 1
        tail = str(Seq(F[fs + (first - fs) // 3 * 3: stop_i + 3]).translate())
    # junction support across all candidate contigs of this species/gene
    sup = {}
    if key:
        a = key[0] - 1
        fj = F[a - 30:a + 30]
        mj_pos = key[2]                        # 1-based last M base before the junction
        mj = M[mj_pos - 30:mj_pos + 30]
        db = "".join(f">{h}\n{s}\n" for h, s in cand.items() if h.split("|")[0] == sp and h.split("|")[3] == gene)
        for name, q in (("female", fj), ("male", mj)):
            hits = blastn(q, db, short=True)
            ok = {x[6] for x in hits if float(x[4]) >= 95 and int(x[5]) >= 56 and int(x[0]) <= 3 and int(x[1]) >= 58}
            per = collections.Counter(h.split("|")[1] for h in ok)
            nF = sum(1 for s_ in per if sex_of[s_] == "female"); nM = sum(1 for s_ in per if sex_of[s_] == "male")
            sup[name] = f"{nF}/{nM}"
            for s_ in sorted(sex_of):
                if s_ in {h.split('|')[1] for h in cand if h.split('|')[0] == sp}:
                    J.write(f"{sp}\t{gene}\t{iso}\t{name}\t{s_}\t{sex_of[s_]}\t{per.get(s_, 0)}\n")
    seg_txt = ";".join(f"{a}-{b}({b - a + 1}bp){'*' if (a, b, c) == key else ''}" for a, b, c in segs)
    P.write(f"{sp}\t{gene}\t{iso}\t{fh}\t{mh}\t{len(fprot) - 1}\t{len(mprot) - 1}\t{seg_txt}\t{cod[-1] if stop_i else ''}\t"
            f"{stop_i + 1 if stop_i is not None else ''}\t{in_spec}\t{ctx}\t{nt4}\t{nxt}\t{nxtc}\t"
            f"{','.join(f'{k}:{v}' for k, v in sorted(collections.Counter(inframe).items()))}\t{tail}\t"
            f"{sup.get('female', '')}\t{sup.get('male', '')}\n")
    fid, mid = fh.split("|")[1] + "_" + fh.split("|")[-1], mh.split("|")[1] + "_" + mh.split("|")[-1]
    TM.write(f">{sp}|{gene}|female|{iso}|{fid} CDS={fs + 1 if fs is not None else 'NA'}-{stop_i + 3 if stop_i is not None else 'NA'} "
             f"female_specific={seg_txt}\n{F}\n")
    TP.write(f">{sp}|{gene}|female|{iso}|{fid}\n{fprot}\n")
    if iso.endswith("1") or iso == "canonical":
        TM.write(f">{sp}|{gene}|male|{mid} CDS={ms + 1}-{me}\n{M}\n")
        TP.write(f">{sp}|{gene}|male|{mid}\n{mprot}\n")
for h in (P, S, J, TM, TP):
    h.close()
for x in os.listdir(tmp):
    os.remove(os.path.join(tmp, x))
os.rmdir(tmp)
