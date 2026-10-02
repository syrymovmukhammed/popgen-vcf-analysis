#!/usr/bin/env bash
# Step 04 — pairwise SNP distance matrix (+ heatmap; + validation against the
# true distances in sim mode).
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3
need_file "${VCF_DIR}/core_snps.vcf.gz"
truth=()
[[ "${MODE}" == "sim" ]] && truth=(--truth "${DATA_DIR}/truth_snp_distance.tsv")
python3 "${PROJECT_DIR}/scripts/snp_distance.py" --vcf "${VCF_DIR}/core_snps.vcf.gz" \
    --outdir "${DIST_DIR}" "${truth[@]}"
