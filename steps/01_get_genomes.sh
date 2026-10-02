#!/usr/bin/env bash
# Step 01 — strain genomes.
#   MODE=sim  : evolve ${N_CLADES} x ${STRAINS_PER_CLADE} genomes from the
#               reference along a known tree (truth labels + truth distances)
#   MODE=real : download the complete chromosomes listed in
#               config/real_genomes.tsv from NCBI (labels = literature phylogroups)
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3
need_file "${REF}"
mkdir -p "${GENOME_DIR}"

if [[ "${MODE}" == "sim" ]]; then
  python3 "${PROJECT_DIR}/scripts/simulate_population.py" --reference "${REF}" --outdir "${DATA_DIR}" \
      --n-clades "${N_CLADES}" --strains-per-clade "${STRAINS_PER_CLADE}" --seed "${SEED}"
else
  need_cmd curl
  printf 'strain\tgroup\n' > "${LABELS}"
  while IFS=$'\t' read -r acc strain group _; do
    [[ "${acc}" == "accession" || -z "${acc}" ]] && continue
    out="${GENOME_DIR}/${strain//[:\/]/_}.fa"
    if [[ ! -s "${out}" ]]; then
      log "downloading ${acc} (${strain}, phylogroup ${group})"
      curl -fsSL --retry 4 --retry-delay 5 --max-time 600 -o "${out}.tmp" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
      { echo ">${strain//[:\/]/_}"; grep -v '>' "${out}.tmp"; } > "${out}"
      rm -f "${out:?}.tmp"
      sleep 0.4   # NCBI asks for <= 3 requests/second without an API key
    fi
    printf '%s\t%s\n' "${strain//[:\/]/_}" "${group}" >> "${LABELS}"
  done < "${REAL_GENOMES_TSV}"
fi
log "$(strains | wc -l) genomes in ${GENOME_DIR}"
