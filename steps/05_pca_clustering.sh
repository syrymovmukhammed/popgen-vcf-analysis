#!/usr/bin/env bash
# Step 05 — PCA, hierarchical + k-means clustering, neighbour-joining tree,
# agreement with known groups; then the run summary (results/summary.md).
source "$(dirname "$0")/../lib/common.sh"
need_cmd python3
need_file "${VCF_DIR}/core_snps.vcf.gz"
python3 "${PROJECT_DIR}/scripts/pca_clustering.py" --vcf "${VCF_DIR}/core_snps.vcf.gz" \
    --labels "${LABELS}" --outdir "${PCA_DIR}" --seed "${SEED}"
python3 "${PROJECT_DIR}/scripts/make_report.py" --project "${PROJECT_DIR}" --out "${RESULTS_DIR}/summary.md"
log "summary -> ${RESULTS_DIR}/summary.md"
