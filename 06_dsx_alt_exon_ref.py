#!/usr/bin/env python
"""Stop codon of the alternative dsx female exon (Tox isoform F3, 'VEMLFHNN...') in the three reference genomes.
The exon is located by tblastn of the Tox F3 peptide against each dsx locus; translation runs from the
start of the hit, in its frame, to the first stop.
Usage: python 06_dsx_alt_exon_ref.py ref_genes.tsv out.tsv
"""
import sys, subprocess, tempfile, os
from Bio.Seq import Seq

PEP = "VEMLFHNNVEMQPFWNYHV"
loci = {"Aedes_aegypti": "AaegL5_1:233156790-233676622",
        "Culex_quinquefasciatus": "NC_051861.1:46625812-46878888",
        "Anopheles_gambiae": "NC_064601.1:47610877-47700920"}
genomes = {l.split("\t")[0]: l.rstrip("\n").split("\t")[4] for l in open(sys.argv[1]) if not l.startswith("species")}
out = open(sys.argv[2], "w")
out.write("species\tgene\tisoform\thit_genomic\tpident\tfemale_specific_peptide\tstop_codon\tstop_genomic\tstop_context_-6_+6\tnt4_after_stop\tcodons_to_next_inframe_stop\tnext_stop_codon\n")
d = tempfile.mkdtemp()
open(f"{d}/q.fa", "w").write(f">q\n{PEP}\n")
for sp, loc in loci.items():
    seq = "".join(subprocess.run(["samtools", "faidx", genomes[sp], loc], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
    open(f"{d}/s.fa", "w").write(f">s\n{seq}\n")
    h = subprocess.run(["tblastn", "-query", f"{d}/q.fa", "-subject", f"{d}/s.fa", "-seg", "no", "-evalue", "0.01",
                        "-outfmt", "6 qstart qend sstart send pident evalue"], capture_output=True, text=True, check=True).stdout.split("\n")[0].split("\t")
    qs, ss, se, pid = int(h[0]), int(h[2]), int(h[3]), h[4]
    chrom, rng = loc.split(":"); lo = int(rng.split("-")[0])
    if ss > se:                                    # minus strand: work on the reverse complement
        seq = str(Seq(seq).reverse_complement()); L = len(seq)
        ss, se = L - ss + 1, L - se + 1
        g = lambda i: lo + (L - 1 - i)            # 0-based rc index -> genomic
    else:
        g = lambda i: lo + i
    i = ss - 1 - 3 * (qs - 1)
    cod = []
    while seq[i:i + 3] not in ("TAA", "TAG", "TGA"):
        cod.append(seq[i:i + 3]); i += 3
    stop = seq[i:i + 3]
    j = i + 3
    while seq[j:j + 3] not in ("TAA", "TAG", "TGA"):
        j += 3
    out.write(f"{sp}\tdsx\tF_alt(5b-type)\t{chrom}:{g(ss - 1)}\t{pid}\t{Seq(''.join(cod)).translate()}*\t{stop}\t{chrom}:{g(i)}\t"
              f"{seq[i - 6:i + 9]}\t{seq[i + 3]}\t{(j - i) // 3}\t{seq[j:j + 3]}\n")
out.close()
for x in os.listdir(d):
    os.remove(os.path.join(d, x))
os.rmdir(d)
