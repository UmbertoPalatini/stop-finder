#!/bin/bash
# BLAST databases of the per-sample de novo Trinity assemblies (Toxorhynchites, Wyeomyia) and
# tblastn of the female/male MOSN, dsx and fru proteins of Aedes, Culex and Anopheles against them.
#SBATCH --job-name=yt_blast
#SBATCH --partition=voss_a
#SBATCH --account=voss_condo_bank
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/03_blast_%j.out
#SBATCH --error=logs/03_blast_%j.err
set -eo pipefail
R=/lustre/fs8/voss_lab/scratch/upalatini/projects/Collabs/Yael/runs/2026-10-29_reviewer_transcripts_stops
T=/lustre/fs8/voss_lab/scratch/upalatini/projects/Collabs/Yael/Trinity
cd $R
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py310
tblastn -version > logs/03_versions.txt; makeblastdb -version >> logs/03_versions.txt
mkdir -p work/blastdb results/03_blast
Q=results/02_ref_stops/ref_FM_protein.fa
for fa in $T/toxo_ambo_brains/trinity_assemblies/*.Trinity.fasta $T/wyeomyia_smithii_brains/trinity_assemblies/*.Trinity.fasta; do
  s=$(basename $fa _trinity.Trinity.fasta)
  [[ -s work/blastdb/$s.nsq || -s work/blastdb/$s.00.nsq ]] || makeblastdb -in $fa -dbtype nucl -parse_seqids -out work/blastdb/$s > /dev/null
  tblastn -query $Q -db work/blastdb/$s -evalue 1e-5 -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 200 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sframe" > results/03_blast/$s.tblastn.tsv
  echo "$s $(wc -l < results/03_blast/$s.tblastn.tsv) hits"
done
echo done $(date)
