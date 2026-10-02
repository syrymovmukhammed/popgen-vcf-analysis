#!/usr/bin/env bash
# Step 03 — VCF QC and core-genome SNP filter (bcftools):
#   * biallelic SNPs only            (-m2 -M2 -v snps)
#   * F_MISSING <= MAX_MISSING       (site present in >= 90% of strains)
#   * still polymorphic after filtering (at least one ALT and one REF call)
# + per-strain missingness / private-SNP stats from `bcftools stats -s -`.
source "$(dirname "$0")/../lib/common.sh"
need_cmd bcftools
in="${VCF_DIR}/all_sites.vcf.gz"
core="${VCF_DIR}/core_snps.vcf.gz"
need_file "${in}"

bcftools view -m2 -M2 -v snps -i "F_MISSING<=${MAX_MISSING}" -Ou "${in}" \
  | bcftools view -i 'COUNT(GT="alt")>0 && COUNT(GT="ref")>0' -Oz -o "${core}"
bcftools index -t "${core}"

bcftools stats -s - "${in}"   > "${VCF_DIR}/all_sites.bcftools_stats.txt"
bcftools stats -s - "${core}" > "${VCF_DIR}/core_snps.bcftools_stats.txt"
# per-strain: number of missing genotypes in the core set
bcftools query -f '[%SAMPLE\t%GT\n]' "${core}" \
  | awk -F'\t' '{n[$1]++; if($2==".") m[$1]++} END{print "strain\tcore_sites\tmissing"; for(s in n) print s"\t"n[s]"\t"(m[s]+0)}' \
  | sort > "${VCF_DIR}/per_strain_missing.tsv"

n_all=$(bcftools view -H "${in}" | wc -l); n_core=$(bcftools view -H "${core}" | wc -l)
log "core SNPs: ${n_core} of ${n_all} variable sites kept (biallelic, missing <= ${MAX_MISSING})"
