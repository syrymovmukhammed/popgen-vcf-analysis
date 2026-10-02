#!/usr/bin/env bash
# Acceptance test used by CI. Fails (exit 1) if an artefact is missing, the
# VCF is malformed, or — in sim mode — if the pipeline does not recover the
# simulated population (distance correlation, clustering ARI).
source "$(dirname "$0")/../lib/common.sh"
fail=0
check() { if eval "$2"; then echo "  ok    $1"; else echo "  FAIL  $1"; fail=1; fi; }

echo "== artefacts"
for f in "${VCF_DIR}/all_sites.vcf.gz" "${VCF_DIR}/core_snps.vcf.gz" "${VCF_DIR}/core_snps.vcf.gz.tbi" \
         "${DIST_DIR}/snp_distance_matrix.tsv" "${DIST_DIR}/snp_distance_heatmap.png" \
         "${PCA_DIR}/clusters.tsv" "${PCA_DIR}/pca_clusters.png" "${PCA_DIR}/nj_tree.nwk" \
         "${PCA_DIR}/nj_tree.png" "${RESULTS_DIR}/summary.md"; do
  check "$(basename "$f")" "[[ -s '$f' ]]"
done

echo "== integrity"
n_strains=$(strains | wc -l)
check "VCF has one column per genome (${n_strains})" \
      "[[ \$(bcftools query -l '${VCF_DIR}/core_snps.vcf.gz' | wc -l) -eq ${n_strains} ]]"
check "core VCF: only biallelic SNPs" \
      "[[ -z \$(bcftools view -H -e 'TYPE=\"snp\" && N_ALT=1' '${VCF_DIR}/core_snps.vcf.gz') ]]"
check "core VCF has >= 100 SNPs" "[[ \$(bcftools view -H '${VCF_DIR}/core_snps.vcf.gz' | wc -l) -ge 100 ]]"
check "distance matrix is ${n_strains}x${n_strains} and symmetric" \
      "python3 -c \"import pandas as p; d=p.read_csv('${DIST_DIR}/snp_distance_matrix.tsv',sep='\t',index_col=0); exit(not (d.shape==(${n_strains},${n_strains}) and (d.values==d.values.T).all()))\""
check "NJ tree contains every strain" \
      "[[ \$(grep -o '[A-Za-z0-9_.-]*:' '${PCA_DIR}/nj_tree.nwk' | grep -c '[A-Za-z]') -ge ${n_strains} ]]"

if [[ "${MODE}" == "sim" ]]; then
  echo "== recovery of the simulated population"
  r=$(python3 -c "import json; print(json.load(open('${DIST_DIR}/distance_stats.json'))['pearson_r_vs_truth'])")
  ari=$(python3 -c "import json; print(json.load(open('${PCA_DIR}/clustering_stats.json'))['ARI_hierarchical'])")
  check "distance vs truth r=${r} >= ${MIN_DIST_R}" "awk 'BEGIN{exit !(${r} >= ${MIN_DIST_R})}'"
  check "clustering ARI=${ari} >= ${MIN_ARI}" "awk 'BEGIN{exit !(${ari} >= ${MIN_ARI})}'"
fi

[[ ${fail} -eq 0 ]] && echo "ALL CHECKS PASSED" || { echo "SOME CHECKS FAILED"; exit 1; }
