# stopfinder

Find the termination codon that ends a female-specific (sex-specifically spliced) isoform, and describe
it: which stop codon it is, whether it lies in the female-specific exon or is caused by a frameshift, its
sequence context, the nucleotide after it, and how far away the next in-frame stop is.

stopfinder works for any species and any gene with a female and a male isoform that share a start codon,
such as *doublesex*, *fruitless* or other sex-specifically spliced genes. It was written for the analysis of
MOSN, *doublesex* and *fruitless* in five mosquito species (see [`paper/`](paper/)), and generalised so that
anyone can run it on their own data.

## What it does

For each female/male isoform pair, stopfinder:

1. defines the **female-specific sequence**: the part of the female mRNA that is absent from the male mRNA;
2. translates the female mRNA from the **start codon shared with the male isoform** to the first in-frame stop;
3. reports, for that stop:
   - the codon (DNA and RNA: UAA, UAG, UGA);
   - where it lies: in the female-specific sequence, across a junction, or in shared sequence downstream
     (when a female exon that is not a multiple of 3 nt shifts the reading frame), and the frame relative to
     the male reading frame;
   - the nucleotide immediately after it (+4) and the sequence context (6 nt on each side);
   - the number of codons to the next in-frame stop, and that codon;
   - the female-specific peptide (from the start of the female-specific sequence to the stop);
4. lists every stop codon in the female-specific sequence, in all three frames;
5. writes 60-nt probes for the female and male splice junctions, which the `junctions` command can search
   in per-sample assemblies to record in which samples each splice form was assembled.

There are two input modes. Both produce the same output table.

| Mode | Input | Female-specific sequence |
|---|---|---|
| `transcripts` | female and male mRNA in FASTA, from any source (de novo assembly, Iso-Seq, NCBI RefSeq mRNAs, ...) | female positions that do not align to the male mRNA |
| `genome` | genome FASTA, GTF or GFF3, and two transcript IDs (plus optional extra exons) | genomic positions in the female transcript that are not in the male transcript |

## Installation

Requires Python ≥ 3.9 and Biopython. No external programs are needed.

```bash
pip install git+https://github.com/UmbertoPalatini/stop-finder.git
```

or from a clone:

```bash
git clone https://github.com/UmbertoPalatini/stop-finder.git
```

```bash
cd stop-finder && pip install .
```

Add `openpyxl` (`pip install ".[xlsx]"`) to also write Excel output with `--xlsx`.
`python -m stopfinder` works in place of the `stopfinder` command.

## Quick start

The `examples/` folder has a small synthetic gene with two female-specific exons: exon A contains a UGA
followed immediately by UAG; exon B is 220 nt long, has no stop and shifts the reading frame.

```bash
cd examples
```

```bash
stopfinder genome --genome genome.fa --annotation genes.gtf --pairs pairs_genome.tsv -o out_genome
```

```bash
stopfinder transcripts --fasta transcripts.fa --pairs pairs_transcripts.tsv -o out_transcripts
```

A single pair needs no table:

```bash
stopfinder transcripts --female female.fa --male male.fa --name dsx -o out
```

## Input

### Transcript mode

`--pairs` is a tab-separated table with a header. Sequence IDs are the first word of the FASTA headers.

| Column | Required | Meaning |
|---|---|---|
| `pair_id` | yes | name of the pair in the output |
| `female` | yes | ID of the female mRNA |
| `male` | yes | ID of the male mRNA |
| `male_start` | no | 1-based position of the first base of the male start codon |
| `female_start` | no | 1-based position of the first base of the shared start codon in the female mRNA |

Both sequences must be in sense (5′→3′) orientation. The male start codon is taken from, in order:
`male_start`; a `CDS=start-end` tag in the male FASTA header; the longest ATG-initiated ORF of the male mRNA.
The female start is the female position aligned to the male start codon, unless `female_start` is given.

Options: `--min-segment` (default 15): unaligned female regions shorter than this are treated as
differences between individuals, not as female-specific sequence.

### Genome mode

`--pairs` is a tab-separated table with a header:

| Column | Required | Meaning |
|---|---|---|
| `pair_id` | yes | name of the pair in the output |
| `female_tx` | yes | female transcript ID in the annotation |
| `male_tx` | yes | male transcript ID in the annotation |
| `extra_female_exons` | no | exons to add to the female transcript, if they are not annotated: `start-end` or `chrom:start-end`, several separated by `;` (`.` for none) |
| `start` | no | genomic coordinate of the first base of the shared start codon (default: the annotated CDS start of `male_tx`) |
| `genome`, `annotation` | no | per-row genome and annotation, to run several species in one table (default: `--genome`, `--annotation`) |

The annotation can be GTF (`transcript_id`) or GFF3 (`Parent=` of exon and CDS features). Transcript IDs
are also matched without the `rna-` and `transcript:` prefixes that NCBI and Ensembl GFF3 files use.
Overlapping exons in the female transcript (for example an added exon that extends an annotated one) are
merged. The genome FASTA must be uncompressed; a samtools-compatible `.fai` index is created next to it if
missing.

If the female isoform is not annotated, use the male transcript ID for both `female_tx` and `male_tx` and
add the female exon with `extra_female_exons`. This works when the female isoform is the male isoform plus
the extra exon (a cassette exon or an exon extension).

### Junction support (`junctions`)

```bash
stopfinder junctions --junctions out/junctions.tsv --samples samples.tsv -o out_junctions
```

`samples.tsv` has columns `sample`, `group` (for example `female`/`male`) and `fasta` (one assembly per
sample; relative paths are read from the folder of `samples.tsv`). For each pair, the female junction is
30 nt before and 30 nt after the start of the female-specific segment that carries the stop; the male
junction is 30 nt on each side of the corresponding male splice site. A sequence contains a junction if
the 60-nt probe matches it on either strand with at most `--max-mismatches` mismatches (default 3,
i.e. ≥ 95% identity) and no gaps. The output gives the number of sequences per sample that contain each
junction, and a summary of how many samples in each group have at least one.

This records whether each splice form was assembled in each sample. It does not count reads.

## Output

| File | Content |
|---|---|
| `summary.tsv` | one row per pair (columns below) |
| `female_specific_stops.tsv` | every stop codon in the female-specific sequence, in all three frames |
| `junctions.tsv` | 60-nt female and male junction probes, input for `stopfinder junctions` |
| `proteins.fa` | female and male proteins |
| `transcripts.fa` | female and male mRNAs, with `CDS=` and `female_specific=` in the headers |
| `summary.xlsx` | with `--xlsx`: the two tables as Excel sheets |

Pairs that cannot be analysed (for example a missing transcript ID) are listed in `summary.tsv` with
`ERROR` in the `note` column; the command then exits with status 1.

### `summary.tsv` columns

| Column | Meaning |
|---|---|
| `female_cds`, `male_cds` | 1-based positions of the start codon and of the last base of the stop codon in each mRNA |
| `start_codon` | codon at the female start position |
| `female_protein_aa`, `male_protein_aa` | protein lengths, without the stop |
| `shared_nterm_aa` | number of identical N-terminal residues of the two proteins |
| `female_specific_segments` | each segment as `start-end(length bp)` in the female mRNA, plus genomic interval in genome mode; `*` marks the segment that carries the stop (the last segment starting before it) |
| `female_specific_bp` | length of that segment |
| `stop_codon`, `stop_codon_rna` | termination codon of the female ORF |
| `stop_pos_female`, `stop_genomic` | position of its first base in the female mRNA, and in the genome (genome mode; transcript orientation, so on the − strand this is the highest coordinate of the codon) |
| `stop_location` | `female-specific` (all three bases), `junction` (one or two bases), `shared, downstream of female-specific sequence`, `shared, upstream of female-specific sequence` (translation ends before the female sequence) or `shared (no female-specific sequence)` |
| `frame_vs_male` | for stops in shared sequence: `male frame`, or `shifted +1`/`+2` (nt offset of the female reading frame from the male one) |
| `nt_plus4` | the nucleotide 3′ of the stop codon |
| `context` | 6 nt upstream, the stop codon, 6 nt downstream |
| `codons_to_next_stop`, `next_stop_codon`, `next_stop_codon_rna` | next in-frame stop after the termination codon (1 = immediately adjacent). This is the longest extension that read-through of the first stop could add |
| `inframe_stops_in_female_specific` | counts of in-frame stops (in the ORF frame) inside the female-specific segments |
| `female_tail_peptide` | translation from the codon containing the first base of the starred segment to the stop |
| `note` | warnings and explanations, for example a frameshift, a non-ATG start, low alignment identity, or how the start codon was chosen |

## Method

The female and male mRNAs are compared directly: in genome mode by genomic coordinate, in transcript mode
by a global pairwise alignment (Biopython `PairwiseAligner`; match +1, mismatch −4, gap open −10, gap extend
−0.1, free end gaps). This scoring aligns near-identical shared exons despite a few SNPs and leaves inserted
or replaced exons unaligned. Aligned islands shorter than 30 nt surrounded by unaligned sequence are treated
as chance matches. Unaligned regions at the 5′ or 3′ end of the female mRNA count as female-specific only
if the male mRNA also has unaligned sequence at that end (an alternative first or last exon); otherwise they
are differences in UTR or contig length.

The female mRNA is then read in codons from the shared start codon to the first stop codon of the chosen
genetic code (`--table`, NCBI numbering, default 1).

## Limits

- **Choosing the isoforms is up to you.** Many genes have several promoters and isoforms. stopfinder
  analyses the pair you give; it does not decide which transcript is "female" and which is "male".
- **Transcript mode cannot see exon boundaries.** When the end of the female exon has the same bases as the
  start of the next exon, the boundary of the female-specific segment is ambiguous by a few nucleotides, and
  the first residue of `female_tail_peptide` can change. The stop codon and its context are not affected
  unless the stop lies within a few nucleotides of the boundary. Genome mode uses annotated exon boundaries.
- In transcript mode, both sequences must be in sense orientation, and the alignment needs memory roughly
  proportional to the product of the two lengths (fine for mRNAs of a few kb; a 8 kb pair takes ~2 s).
- Read-through is not modelled; `codons_to_next_stop` gives the information needed to reason about it.

## Tests

```bash
pip install ".[test]"
```

```bash
pytest
```

The tests use synthetic genes with known answers (cassette, frameshifting, terminal and mutually exclusive
female exons; both strands; GTF and GFF3; SNPs and UTR differences between female and male sequences), and
check that genome and transcript mode agree.

## Citation

If you use stopfinder, please cite it as given in [`CITATION.cff`](CITATION.cff) (GitHub shows this as
"Cite this repository").

## License

MIT, see [`LICENSE`](LICENSE).
