#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# All parameters in one place; override from the environment, e.g.
#     MODE=real THREADS=8 bash run_all.sh
# ---------------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${THREADS:=2}"
: "${SEED:=2024}"

# MODE=sim  (default, CI): genomes are simulated along a known tree with four
#                          phylogroup-like clades, so clustering can be scored.
# MODE=real               : complete genomes of real E. coli strains listed in
#                          config/real_genomes.tsv are downloaded from NCBI.
: "${MODE:=sim}"

# Reference: E. coli K-12 MG1655 (NC_000913.3). Sim mode uses the first
# REGION_END bp; real mode needs the whole chromosome.
: "${REF_ACC:=NC_000913.3}"
if [[ "${MODE}" == "real" ]]; then : "${REGION_END:=4641652}"; else : "${REGION_END:=200000}"; fi
: "${REF_FASTA:=}"

# Simulation design (sim mode)
: "${N_CLADES:=4}"
: "${STRAINS_PER_CLADE:=6}"
: "${REAL_GENOMES_TSV:=${PROJECT_DIR}/config/real_genomes.tsv}"

# Genome -> reference alignment and site filters
: "${MM2_PRESET:=asm20}"       # minimap2 preset for ~1-3% divergence between E. coli phylogroups
: "${MIN_ALN_LEN:=1000}"       # ignore alignment blocks shorter than this
: "${MAX_MISSING:=0.10}"       # drop sites missing in >10% of strains (core-genome SNPs)

# CI acceptance (sim mode): clustering must recover the simulated clades
: "${MIN_ARI:=0.90}"
: "${MIN_DIST_R:=0.98}"

DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results"
REF="${DATA_DIR}/reference.fa"
GENOME_DIR="${DATA_DIR}/genomes"
LABELS="${DATA_DIR}/strain_labels.tsv"     # strain -> group (truth clade or expected phylogroup)
ALN_DIR="${RESULTS_DIR}/02_alignments"
VCF_DIR="${RESULTS_DIR}/03_vcf"
DIST_DIR="${RESULTS_DIR}/04_distance"
PCA_DIR="${RESULTS_DIR}/05_pca_clustering"

export THREADS SEED MODE
