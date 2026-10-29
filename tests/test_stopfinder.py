import csv
import os
import random

import pytest

from stopfinder import cli, transcripts
from stopfinder.core import analyse
from stopfinder.seqio import IndexedFasta, revcomp

from synth import build_genome, coding, gene_model, mrnas, noise, stopless

KEYS = ["start_codon", "female_protein_aa", "male_protein_aa", "female_specific_bp", "stop_codon",
        "stop_pos_female", "stop_location", "frame_vs_male", "nt_plus4", "context", "codons_to_next_stop",
        "next_stop_codon", "inframe_stops_in_female_specific", "female_tail_peptide"]


@pytest.fixture(scope="module")
def model():
    return gene_model()


@pytest.fixture(scope="module")
def genome_files(model, tmp_path_factory):
    d = tmp_path_factory.mktemp("genome")
    fasta, gtf, coords, L = build_genome(model)
    (d / "genome.fa").write_text(fasta)
    (d / "genes.gtf").write_text(gtf)
    (d / "pairs.tsv").write_text(
        "pair_id\tfemale_tx\tmale_tx\n"
        "A_plus\tFA\tM\nB_plus\tFB\tM\nA_minus\tFA.minus\tM.minus\nB_minus\tFB.minus\tM.minus\n")
    return d, coords, L


def read_rows(path):
    with open(path) as fh:
        return {r["pair_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def run_genome(d, out):
    rc = cli.main(["genome", "--genome", str(d / "genome.fa"), "--annotation", str(d / "genes.gtf"),
                   "--pairs", str(d / "pairs.tsv"), "-o", str(out)])
    assert rc == 0
    return read_rows(out / "summary.tsv")


def test_double_stop_in_female_exon(model, genome_files, tmp_path):
    d, coords, L = genome_files
    rows = run_genome(d, tmp_path)
    r = rows["A_plus"]
    up = len(model["e1"]) + len(model["e2"])
    assert r["stop_codon"] == "TGA" and r["stop_codon_rna"] == "UGA"
    assert int(r["stop_pos_female"]) == up + 60 + 1
    assert r["stop_location"] == "female-specific"
    assert r["nt_plus4"] == "T"
    assert r["codons_to_next_stop"] == "1" and r["next_stop_codon"] == "TAG"
    assert int(r["female_protein_aa"]) == 41 + 50 + 20
    assert int(r["male_protein_aa"]) == 41 + 50 + 60 + 30
    assert int(r["female_specific_bp"]) == 96
    assert r["inframe_stops_in_female_specific"] == "TAG:1,TGA:1"
    s, e = coords["efa"]
    assert r["female_specific_segments"] == f"{up + 1}-{up + 96}(96bp) chr1:{s}-{e}*"
    assert r["stop_genomic"] == f"chr1:{s + 60}(+)"


def test_frameshift_exon_without_stop(model, genome_files, tmp_path):
    d, _, _ = genome_files
    r = run_genome(d, tmp_path)["B_plus"]
    male, fa, fb = mrnas(model)
    fs = model["utr5"]
    k = next(j for j in range(fs, len(fb) - 2, 3) if fb[j:j + 3] in {"TAA", "TAG", "TGA"})
    seg_end = len(model["e1"]) + len(model["e2"]) + 220
    assert k >= seg_end                       # the stop lies after the female exon
    assert int(r["stop_pos_female"]) == k + 1
    assert r["stop_location"] == "shared, downstream of female-specific sequence"
    assert r["frame_vs_male"].startswith("shifted")
    assert r["inframe_stops_in_female_specific"] == "none"
    assert "not a multiple of 3" in r["note"] and f"{k - seg_end} nt into" in r["note"]


def test_minus_strand_matches_plus_strand(genome_files, tmp_path):
    d, coords, L = genome_files
    rows = run_genome(d, tmp_path)
    for case in ("A", "B"):
        p, m = rows[f"{case}_plus"], rows[f"{case}_minus"]
        for k in KEYS:
            assert p[k] == m[k], k
    s, e = coords["efa"]
    assert rows["A_minus"]["stop_genomic"] == f"chr2:{L - (s + 60) + 1}(-)"


def test_transcript_mode_matches_genome_mode(model, genome_files, tmp_path):
    d, _, _ = genome_files
    g = run_genome(d, tmp_path / "g")
    male, fa, fb = mrnas(model)
    (tmp_path / "seqs.fa").write_text(f">M\n{male}\n>FA\n{fa}\n>FB\n{fb}\n")
    (tmp_path / "pairs.tsv").write_text("pair_id\tfemale\tmale\tmale_start\nA_plus\tFA\tM\t41\nB_plus\tFB\tM\t41\n")
    assert cli.main(["transcripts", "--pairs", str(tmp_path / "pairs.tsv"), "--fasta", str(tmp_path / "seqs.fa"),
                     "-o", str(tmp_path / "t")]) == 0
    t = read_rows(tmp_path / "t" / "summary.tsv")
    # Without exon coordinates the female-exon boundary can slide by a base or two when the end of the
    # exon repeats the start of the next one, which can change the first residue of the tail peptide.
    for pid in ("A_plus", "B_plus"):
        for k in KEYS:
            if k == "female_tail_peptide":
                assert t[pid][k].endswith(g[pid][k][1:]), pid
            else:
                assert t[pid][k] == g[pid][k], (pid, k)


def test_transcript_mode_infers_start_and_tolerates_snps_and_utrs(model):
    rng = random.Random(7)
    male, fa, _ = mrnas(model)
    # female contig from another individual: SNPs in shared sequence, shorter 5' end, longer 3' end
    f = list(fa[20:] + noise(150, rng))
    up = len(model["e1"]) + len(model["e2"]) - 20
    for i in (60, 200, up + 150, up + 250):
        f[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[f[i]]
    F = "".join(f)
    p = transcripts.build_pair("x", "F", F, "M", male)
    row = analyse(p)[0]
    assert len(p.segments) == 1                     # the UTR difference is not called female-specific
    a, b = p.segments[0]
    assert b - a == 96 and abs(a - up) <= 2
    assert p.female_start == model["utr5"] - 20
    assert row["stop_codon"] == "TGA" and row["stop_location"] == "female-specific"
    assert row["codons_to_next_stop"] == 1


def test_alternative_terminal_female_exon():
    """dsx-like: the female isoform ends in its own last exon; the male continues with male exons."""
    rng = random.Random(3)
    shared = stopless(30, rng) + "ATG" + coding(80, rng)
    f_last = coding(15, rng) + "TGA" + coding(3, rng) + "TAG" + noise(200, rng)
    m_last = coding(70, rng) + "TAA" + noise(150, rng)
    p = transcripts.build_pair("t", "F", shared + f_last, "M", shared + m_last)
    row = analyse(p)[0]
    # bases shared by chance at the start of the two last exons stay aligned, so allow a few nt
    assert len(p.segments) == 1
    a, b = p.segments[0]
    assert 0 <= a - len(shared) <= 6 and b == len(shared + f_last)
    assert row["stop_location"] == "female-specific" and row["stop_codon"] == "TGA"
    assert row["codons_to_next_stop"] == 4 and row["next_stop_codon"] == "TAG"


def test_cds_tag_in_header(model, tmp_path):
    male, fa, _ = mrnas(model)
    (tmp_path / "f.fa").write_text(f">FA\n{fa}\n")
    (tmp_path / "m.fa").write_text(f">M CDS=41-583\n{male}\n")
    assert cli.main(["transcripts", "--female", str(tmp_path / "f.fa"), "--male", str(tmp_path / "m.fa"),
                     "-o", str(tmp_path / "o")]) == 0
    r = read_rows(tmp_path / "o" / "summary.tsv")["pair1"]
    assert r["female_cds"].startswith("41-") and "longest ATG" not in r["note"]


def test_junction_support(model, tmp_path):
    male, fa, _ = mrnas(model)
    p = transcripts.build_pair("A", "FA", fa, "M", male)
    probes = analyse(p)[2]
    assert {x["junction"] for x in probes} == {"female", "male"}
    (tmp_path / "jn.tsv").write_text("pair_id\tjunction\tsequence\n" +
                                     "".join(f"A\t{x['junction']}\t{x['sequence']}\n" for x in probes))
    rng = random.Random(5)
    (tmp_path / "f1.fa").write_text(f">c1\n{revcomp(fa)}\n>c2\n{noise(500, rng)}\n")   # reverse strand
    (tmp_path / "m1.fa").write_text(f">c1\n{male}\n")
    (tmp_path / "samples.tsv").write_text("sample\tgroup\tfasta\nf1\tfemale\tf1.fa\nm1\tmale\tm1.fa\n")
    assert cli.main(["junctions", "--junctions", str(tmp_path / "jn.tsv"), "--samples",
                     str(tmp_path / "samples.tsv"), "-o", str(tmp_path / "o")]) == 0
    s = {(r["pair_id"], r["junction"]): r["samples_with_junction"]
         for r in csv.DictReader(open(tmp_path / "o" / "junction_support_summary.tsv"), delimiter="\t")}
    assert s[("A", "female")] == "female: 1/1; male: 0/1"
    assert s[("A", "male")] == "female: 0/1; male: 1/1"


def test_indexed_fasta(tmp_path):
    rng = random.Random(9)
    a, b = noise(1234, rng), noise(77, rng)
    (tmp_path / "g.fa").write_text(">a desc\n" + "\n".join(a[i:i + 50] for i in range(0, len(a), 50)) +
                                   "\n>b\n" + b + "\n")
    fa = IndexedFasta(str(tmp_path / "g.fa"))
    assert fa.fetch("a", 1, 1234) == a
    assert fa.fetch("a", 49, 103) == a[48:103]
    assert fa.fetch("b", 70, 77) == b[69:]
    assert os.path.exists(tmp_path / "g.fa.fai")
    assert IndexedFasta(str(tmp_path / "g.fa")).fetch("a", 100, 160) == a[99:160]


def test_gff3_and_extra_exons(model, genome_files, tmp_path):
    """GFF3 input with NCBI-style 'rna-' IDs, and a female exon added by coordinates."""
    d, coords, _ = genome_files
    gff = []
    for line in (d / "genes.gtf").read_text().splitlines():
        f = line.split("\t")
        tid = f[8].split('transcript_id "')[1].split('"')[0]
        if f[0] == "chr1" and tid in ("M", "FA"):
            f[8] = f"Parent=rna-{tid}"
            gff.append("\t".join(f))
    (tmp_path / "g.gff3").write_text("##gff-version 3\n" + "\n".join(gff) + "\n")
    s, e = coords["efa"]
    (tmp_path / "pairs.tsv").write_text(
        "pair_id\tfemale_tx\tmale_tx\textra_female_exons\nA\tFA\tM\t.\nA_added\tM\tM\t" + f"chr1:{s}-{e}\n")
    assert cli.main(["genome", "--genome", str(d / "genome.fa"), "--annotation", str(tmp_path / "g.gff3"),
                     "--pairs", str(tmp_path / "pairs.tsv"), "-o", str(tmp_path / "o"), "--xlsx"]) == 0
    rows = read_rows(tmp_path / "o" / "summary.tsv")
    for k in KEYS:
        assert rows["A"][k] == rows["A_added"][k], k
    assert rows["A"]["stop_codon"] == "TGA"
    assert (tmp_path / "o" / "summary.xlsx").exists()


def test_failed_pair_is_reported(genome_files, tmp_path):
    d, _, _ = genome_files
    (tmp_path / "pairs.tsv").write_text("pair_id\tfemale_tx\tmale_tx\nbad\tNOPE\tM\nA\tFA\tM\n")
    rc = cli.main(["genome", "--genome", str(d / "genome.fa"), "--annotation", str(d / "genes.gtf"),
                   "--pairs", str(tmp_path / "pairs.tsv"), "-o", str(tmp_path / "o")])
    assert rc == 1
    rows = read_rows(tmp_path / "o" / "summary.tsv")
    assert rows["bad"]["note"].startswith("ERROR") and rows["A"]["stop_codon"] == "TGA"


@pytest.mark.parametrize("seed", range(60))
def test_random_gene_structures(seed):
    from synth import random_pair
    kind, F, M, (a, b) = random_pair(seed)
    p = transcripts.build_pair("r", "F", F, "M", M)
    row = analyse(p)[0]
    assert row["stop_location"] == "female-specific", (kind, p.segments, (a, b))
    assert len(p.segments) == 1
    assert abs(p.segments[0][0] - a) <= 6 and abs(p.segments[0][1] - b) <= 6


def test_single_pair_with_identical_contig_names(model, tmp_path):
    """Trinity names restart in every assembly, so a female and a male contig can share a name."""
    male, fa, _ = mrnas(model)
    (tmp_path / "f.fa").write_text(f">TRINITY_DN1_c0_g1_i1\n{fa}\n")
    (tmp_path / "m.fa").write_text(f">TRINITY_DN1_c0_g1_i1\n{male}\n")
    assert cli.main(["transcripts", "--female", str(tmp_path / "f.fa"), "--male", str(tmp_path / "m.fa"),
                     "-o", str(tmp_path / "o")]) == 0
    assert read_rows(tmp_path / "o" / "summary.tsv")["pair1"]["stop_codon"] == "TGA"


def test_uneven_fasta_lines_are_rejected(tmp_path):
    (tmp_path / "g.fa").write_text(">c\n" + "A" * 60 + "\n" + "C" * 40 + "\n" + "G" * 60 + "\n")
    with pytest.raises(ValueError, match="different lengths"):
        IndexedFasta(str(tmp_path / "g.fa"))
    (tmp_path / "ok.fa").write_text(">c\n" + "A" * 60 + "\n" + "C" * 40 + "\n\n>d\nACGT\n")
    assert IndexedFasta(str(tmp_path / "ok.fa")).fetch("c", 55, 65) == "AAAAAACCCCC"


def test_table_row_with_extra_column(tmp_path):
    (tmp_path / "p.tsv").write_text("pair_id\tfemale\tmale\nA\tF\tM\textra\n")
    with pytest.raises(SystemExit, match="more columns than the header"):
        cli.read_table(str(tmp_path / "p.tsv"), ["pair_id"])


def test_female_specific_5utr_does_not_carry_the_stop(model, genome_files, tmp_path):
    """A female-only 5' UTR stretch must not be taken as the segment that carries the stop."""
    d, coords, _ = genome_files
    s, e = coords["e1"]
    gtf = (d / "genes.gtf").read_text()
    (tmp_path / "g.gtf").write_text(gtf + f'chr1\ttest\texon\t{s - 300}\t{s - 201}\t.\t+\t.\tgene_id "G"; '
                                          f'transcript_id "FA";\n')
    (tmp_path / "p.tsv").write_text("pair_id\tfemale_tx\tmale_tx\nA\tFA\tM\n")
    assert cli.main(["genome", "--genome", str(d / "genome.fa"), "--annotation", str(tmp_path / "g.gtf"),
                     "--pairs", str(tmp_path / "p.tsv"), "-o", str(tmp_path / "o")]) == 0
    r = read_rows(tmp_path / "o" / "summary.tsv")["A"]
    assert r["stop_location"] == "female-specific" and r["female_specific_bp"] == "96"
    assert r["female_specific_segments"].startswith("1-100(100bp)")


def test_inframe_count_uses_only_the_segment_with_the_stop(model, genome_files, tmp_path):
    """Female-only sequence in the 5' UTR or downstream of the stop must not enter the in-frame stop count."""
    d, coords, _ = genome_files
    s1 = coords["e1"][0]
    e4 = coords["e4"][1]
    gtf = (d / "genes.gtf").read_text()
    extra = "".join(f'chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\tgene_id "G"; transcript_id "FA";\n'
                    for a, b in ((s1 - 300, s1 - 201), (e4 + 250, e4 + 549)))
    (tmp_path / "g.gtf").write_text(gtf + extra)
    (tmp_path / "p.tsv").write_text("pair_id\tfemale_tx\tmale_tx\nA\tFA\tM\n")
    assert cli.main(["genome", "--genome", str(d / "genome.fa"), "--annotation", str(tmp_path / "g.gtf"),
                     "--pairs", str(tmp_path / "p.tsv"), "-o", str(tmp_path / "o")]) == 0
    r = read_rows(tmp_path / "o" / "summary.tsv")["A"]
    assert r["female_specific_segments"].count("bp") == 3
    assert r["inframe_stops_in_female_specific"] == "TAG:1,TGA:1"
    with open(tmp_path / "o" / "female_specific_stops.tsv") as fh:
        segs = {row["segment"] for row in csv.DictReader(fh, delimiter="\t")}
    assert len(segs) >= 2          # stops of the other segments are still listed in the per-stop table
