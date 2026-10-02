#!/usr/bin/env bash
# Step 02 — align every genome to the reference with minimap2 (assembly
# preset, cs tag on) and convert the alignments into one multi-sample
# haploid VCF (scripts/paf_to_vcf.py). bcftools then sorts, compresses,
# indexes and validates the file.
source "$(dirname "$0")/../lib/common.sh"
need_cmd minimap2 bcftools python3
need_file "${REF}"
mkdir -p "${ALN_DIR}" "${VCF_DIR}"

pafs=()
for s in $(strains); do
  minimap2 -c --cs -x "${MM2_PRESET}" -t "${THREADS}" --secondary=no \
      "${REF}" "${GENOME_DIR}/${s}.fa" > "${ALN_DIR}/${s}.paf" 2> /dev/null
  pafs+=("${ALN_DIR}/${s}.paf")
done

python3 "${PROJECT_DIR}/scripts/paf_to_vcf.py" --reference "${REF}" --paf "${pafs[@]}" \
    --min-aln-len "${MIN_ALN_LEN}" --out "${VCF_DIR}/all_sites.vcf" 2> "${VCF_DIR}/paf_to_vcf.log"
tail -n 3 "${VCF_DIR}/paf_to_vcf.log" >&2
bcftools sort -Oz -o "${VCF_DIR}/all_sites.vcf.gz" "${VCF_DIR}/all_sites.vcf" 2> /dev/null
bcftools index -t "${VCF_DIR}/all_sites.vcf.gz"
rm -f "${VCF_DIR:?}/all_sites.vcf"
log "multi-sample VCF: $(bcftools query -l "${VCF_DIR}/all_sites.vcf.gz" | wc -l) strains, $(bcftools view -H "${VCF_DIR}/all_sites.vcf.gz" | wc -l) variable sites"
