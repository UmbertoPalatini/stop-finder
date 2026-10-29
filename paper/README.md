# Original analysis for the paper (archived)

These are the scripts as they were run for the reviewer response on MOSN, *doublesex* and *fruitless* in
five mosquito species (run `2026-10-29_reviewer_transcripts_stops`). They are kept unchanged so the
published results can be traced to the exact code. For new analyses, use the general `stopfinder` tool
described in the [main README](../README.md); it implements the same stop-codon logic as steps 02 and 05.

## Rerunning the five species with stopfinder

`stopfinder_rerun/` has the same pairs in stopfinder's format. On the HPC, from this folder, with the
original run's `results/` available:

```bash
stopfinder genome --pairs stopfinder_rerun/ref_pairs_genome.tsv -o rerun/genome
```

```bash
stopfinder transcripts --pairs stopfinder_rerun/toxwye_pairs_transcripts.tsv --fasta results/04_classify/candidate_contigs_oriented.fa -o rerun/transcripts
```

```bash
python stopfinder_rerun/compare_with_original.py results rerun/genome rerun/transcripts
```

The last command lists every field (stop codon, +4 nucleotide, context, distance to and identity of the
next stop, protein lengths, in-frame stops, whether the stop is in female-specific sequence) that differs
from the original tables. Two differences are expected, both in `inframe_stops_in_female_specific`:
stopfinder counts in-frame stops only in the female-specific segment that carries the termination codon,
while the original scripts summed all female-specific segments. For *Culex* fru this leaves out a
female-specific 5′ exon upstream of the start codon (16/2/5 instead of 23/5/6); for *Wyeomyia* fru it leaves
out a 1,679 bp female-only 3′ end of the contig, downstream of the stop (12/3/10 instead of 26/4/13). The
stopfinder values are the correct ones. Step 06 (the alternative *doublesex* exon located by a peptide search) and the
contig search (steps 03 and 04) are specific to this study and are not part of stopfinder.

---

## Termination codons in female-specific exons of MOSN, doublesex and fruitless

Code for the analysis of termination codons in the female-specific (male-skipped) exons of
MOSN (AAEL011211 and orthologs), *doublesex* and *fruitless* in five mosquito species:
*Aedes aegypti*, *Culex quinquefasciatus*, *Anopheles gambiae* (reference genomes) and
*Toxorhynchites amboinensis*, *Wyeomyia smithii* (de novo Trinity transcriptome assemblies).

For each gene and species, the scripts take a female and a male transcript and do the following:
- define the female-specific sequence (present in the female transcript, absent from the male one);
- translate the female transcript from the start codon it shares with the male isoform to the first in-frame stop;
- report the termination codon (UAA, UAG or UGA), whether it lies in the female-specific sequence, the
  nucleotide after it (+4), its sequence context, the distance to the next in-frame stop, and every stop
  codon in the female-specific sequence.

### Requirements
- Python 3.10 with Biopython 1.84 and openpyxl 3.1 (openpyxl only for step 07)
- samtools 1.6 (`faidx`), BLAST+ 2.14.0 (`tblastn`, `blastn`, `makeblastdb`, `blastdbcmd`)
- SLURM for step 03 (the script is a plain bash loop and also runs without SLURM)

### Inputs
- Reference genomes (FASTA, indexed with `samtools faidx`) and GTF annotations that include the
  female-specific transcripts: *Ae. aegypti* AaegL5 (VectorBase release 58), *Cx. quinquefasciatus*
  VPISU_Cqui_1.0 (GCF_015732765.1), *An. gambiae* idAnoGambNW_F1_1 (GCF_943734735.2).
- Per-sample Trinity assemblies (Trinity v2.1.1) of male and female brain RNA-seq,
  BioProject PRJNA612100 (3 males and 3 females per species).

The paths in `config/` are absolute paths on the Rockefeller HPC. Edit them before running elsewhere.

### Config files (`config/`)
| File | Content |
|---|---|
| `ref_genes.tsv` | species, gene, gene ID, GTF, genome FASTA for the three reference species |
| `ref_pairs.tsv` | female and male transcript per gene; `add_exons` adds a female-specific exon not present in the GTF (coordinates `start-end`) |
| `samples.tsv` | species, SRA run, sex for the 12 Toxorhynchites and Wyeomyia assemblies |
| `toxwye_pairs.tsv` | representative female and male contig per gene and female isoform (contig names as written by step 04) |

### Steps
| Step | Script | What it does | Main output |
|---|---|---|---|
| 01 | `01_ref_isoforms.py` | Extracts exon structure, mRNA, CDS and protein of every transcript of the three genes from each reference genome and GTF. | `ref_transcripts.tsv`, `ref_exons.tsv`, FASTA files |
| 02 | `02_ref_stops.py` | Reference species. Builds the female transcript (overlapping exons merged), takes the annotated start codon of the male transcript, translates the female transcript to the first stop and reports the termination codon and its context. Female-specific sequence = genomic positions in the female transcript that are absent from the male transcript. | `ref_pairs.tsv`, `ref_female_specific_stops.tsv`, `ref_FM_protein.fa` |
| 03 | `03_blast_assemblies.sbatch` | Builds BLAST databases of the Trinity assemblies and searches them with the female and male proteins from step 02 (tblastn, E ≤ 1e-5). | `<sample>.tblastn.tsv` |
| 04 | `04_classify_contigs.py` | Collects candidate contigs for each gene (≥ 80 bits), orients them, and makes a first female/male call at the sex-specific junction. Used to choose the representative contigs in `toxwye_pairs.tsv`. | `contig_classes.tsv`, `candidate_contigs_oriented.fa` |
| 05 | `05_toxwye_pairs.py` | Toxorhynchites and Wyeomyia. For each female/male contig pair: female-specific sequence = insertion in the female contig relative to the male contig (blastn); the male start codon is located on the female contig by its N-terminal peptide; translation to the first stop and the same stop report as step 02. Also searches the 60-nt female and male junctions in all candidate contigs of the 12 assemblies, to record in which samples each splice form was assembled. | `toxwye_pairs.tsv`, `toxwye_female_specific_stops.tsv`, `toxwye_junction_support.tsv`, `toxwye_transcripts.fa`, `toxwye_proteins.fa` |
| 06 | `06_dsx_alt_exon_ref.py` | Finds the alternative *doublesex* female exon (found in Toxorhynchites transcripts) in the three reference genomes by tblastn of its peptide, and reports its termination codon. | `ref_dsx_alt_exon_stops.tsv` |
| 07 | `07_summary.py` | Combines steps 02, 05 and 06 into one table. | `stop_codon_summary.tsv`, `stop_codon_summary.xlsx` |

### Running
```bash
python 01_ref_isoforms.py config/ref_genes.tsv results/01_ref
python 02_ref_stops.py config/ref_pairs.tsv config/ref_genes.tsv results/01_ref/ref_exons.tsv results/02_ref_stops
sbatch 03_blast_assemblies.sbatch            # edit R and T (run and assembly folders) first
python 04_classify_contigs.py config/samples.tsv results/03_blast work/blastdb results/02_ref_stops/ref_FM_protein.fa results/04_classify
python 05_toxwye_pairs.py config/toxwye_pairs.tsv results/04_classify/candidate_contigs_oriented.fa config/samples.tsv results/05_toxwye
python 06_dsx_alt_exon_ref.py config/ref_genes.tsv results/06_dsx_alt/ref_dsx_alt_exon_stops.tsv
python 07_summary.py results
```
Steps 01, 02 and 06 need `samtools` on the PATH; steps 04–06 need BLAST+.

### Definitions
- **Termination codon:** the first in-frame stop codon when the female transcript is translated from the
  start codon shared with the male isoform.
- **Stop in female-specific sequence:** all three bases of the stop lie in the female-specific sequence. When a
  female exon is not a multiple of three nucleotides long, it can shift the reading frame, so that the stop lies
  in the next shared exon (*Toxorhynchites* MOSN).
- **+4 nucleotide:** the base immediately 3′ of the stop codon.
- **Codons to next in-frame stop:** the number of codons from the first stop to the next stop in the same frame
  (the maximum extension if the first stop were read through).
- **Junction support:** the number of assembled contigs per sample containing the female or male splice junction
  (60 nt, ≥ 95% identity). This records whether each splice form was detected; it does not count reads.
