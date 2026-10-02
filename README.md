# popgen-vcf-analysis — population genomics from a multi-sample VCF

[![CI](https://github.com/syrymovmukhammed/popgen-vcf-analysis/actions/workflows/main.yml/badge.svg)](https://github.com/syrymovmukhammed/popgen-vcf-analysis/actions/workflows/main.yml)
![minimap2](https://img.shields.io/badge/minimap2-asm20-blue.svg)
![python](https://img.shields.io/badge/python-3.12-blue.svg)

From a set of bacterial genomes to population structure:
**genome → reference alignment → multi-sample SNP VCF → core-SNP filter →
SNP-distance matrix → PCA, clustering and a neighbour-joining tree.**
The same code runs on a simulated population with a known answer (CI) or on
15 real *E. coli* genomes spanning the main phylogroups.

## Two modes

| | `MODE=sim` (default, CI) | `MODE=real` |
|---|---|---|
| genomes | 24 strains evolved from the MG1655 sequence along a known tree: 4 clades named after phylogroups A, B1, B2, D, two subclades each | 15 complete chromosomes from NCBI, listed in `config/real_genomes.tsv` (K-12, B, IAI1, 55989, CFT073, UTI89, 536, UMN026, IAI39, O157:H7 Sakai/EDL933 …) |
| reference | first 200 kb of NC_000913.3 | full NC_000913.3 |
| what is checked | distances vs exact simulated distances; clusters vs true clades | clusters vs literature phylogroups |
| runtime | ~10 s | ~10 min |

The simulated population deliberately includes what makes real data awkward:
small indels on every branch, accessory-gene loss (1–4 kb deletions → missing
data) and one **recombination** event that imports a 6 kb segment from another
clade (homoplasy).

## Pipeline

| Step | Script | Tool | Output |
|---|---|---|---|
| 00 | `steps/00_get_reference.sh` | NCBI E-utilities | `data/reference.fa` |
| 01 | `steps/01_get_genomes.sh` | `scripts/simulate_population.py` / NCBI | `data/genomes/*.fa`, `strain_labels.tsv` |
| 02 | `steps/02_genomes_to_vcf.sh` | `minimap2 -c --cs -x asm20` + `scripts/paf_to_vcf.py` + bcftools | `all_sites.vcf.gz` (haploid, one column per strain) |
| 03 | `steps/03_vcf_qc_filter.sh` | bcftools view / stats | `core_snps.vcf.gz`: biallelic SNPs present in ≥ 90 % of strains |
| 04 | `steps/04_snp_distance.sh` | `scripts/snp_distance.py` | distance matrix (snp-dists format), heatmap, truth comparison |
| 05 | `steps/05_pca_clustering.sh` | `scripts/pca_clustering.py` (scikit-learn, SciPy, Biopython) | PCA, UPGMA + k-means clusters, NJ tree (Newick + PNG), ARI/NMI, `results/summary.md` |

How genomes become a VCF (`paf_to_vcf.py`): minimap2's `cs` tag spells out every
difference of an alignment block (`:120` match, `*ag` substitution, `-tt` deletion,
`+ca` insertion). Walking it along the reference gives, per strain, which positions are
covered and which carry a SNP. Positions covered twice (repeats) become missing rather
than fake SNPs; uncovered positions (deletions, gene loss) become `.`. The approach is the
one behind snippy-core / Parsnp core-genome SNP tables, and it is verified to give identical
calls for a genome supplied reverse-complemented and rotated.

Clustering does not assume the number of groups: *k* is picked by the silhouette score
of the UPGMA clustering, then k-means on the leading principal components is run as an
independent second method.

## Running

```bash
sudo apt-get install -y minimap2 samtools bcftools tabix
pip install -r requirements.txt
bash run_all.sh               # MODE=sim
bash tests/check_outputs.sh   # CI acceptance test
MODE=real THREADS=8 bash run_all.sh   # 15 real genomes (needs internet, ~80 MB)
```

The real run can also be started from GitHub: **Actions → CI → Run workflow → mode: real**.

## Results (sim mode, reference local run, offline fallback reference)

- 24 strains, 4,771 core SNPs; measured vs true pairwise distance: **Pearson r = 0.99997**
  (mean absolute error 26 SNPs, from sites lost to strain-specific deletions);
- silhouette picks **k = 4**; ARI vs true clades = **1.0** for both UPGMA and k-means;
- PC1 + PC2 explain 78 % of the variance; the NJ tree recovers every clade and subclade.

Plots: `results/04_distance/snp_distance_heatmap.png`, `distance_vs_truth.png`,
`results/05_pca_clustering/pca_clusters.png`, `nj_tree.png` (also uploaded by CI as an artifact).

## Full-size data

- reference: [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3) (`MODE=real` sets `REGION_END=4641652`);
- genomes: accessions in `config/real_genomes.tsv`, downloaded with
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<ACCESSION>&rettype=fasta`;
  add rows to analyse more strains (bulk alternative: NCBI `datasets download genome taxon "Escherichia coli" --assembly-level complete`).
- any other multi-sample haploid VCF (e.g. from snippy-core) can go straight into steps 03–05:
  copy it to `results/03_vcf/all_sites.vcf.gz` and run `bash steps/03_vcf_qc_filter.sh` onwards.

## Layout

```
popgen-vcf-analysis/
├── .github/workflows/main.yml       # CI: sim run + acceptance test on every push to main
├── config.sh, lib/common.sh, run_all.sh
├── config/real_genomes.tsv       # accession, strain, expected phylogroup
├── steps/00…05_*.sh
├── scripts/
│   ├── simulate_population.py    # genomes on a known tree (+ truth distances)
│   ├── paf_to_vcf.py             # minimap2 cs-tag -> multi-sample haploid VCF
│   ├── vcf_matrix.py             # VCF -> genotype matrix, SNP distance
│   ├── snp_distance.py           # matrix, heatmap, truth comparison
│   ├── pca_clustering.py         # PCA, UPGMA/k-means, silhouette, ARI, NJ tree
│   ├── wgs_sim.py                # FASTA helpers + synthetic fallback reference
│   └── make_report.py
├── tests/check_outputs.sh
└── requirements.txt / environment.yml
```

## Related repositories

Part of a three-repository set on bacterial WGS:

- [wgs-variant-pipeline](https://github.com/syrymovmukhammed/wgs-variant-pipeline) — reads → BAM → VCF, bcftools + GATK, benchmarked against a truth set
- [variant-ml-filter](https://github.com/syrymovmukhammed/variant-ml-filter) — Random Forest / Logistic Regression filter for false-positive calls
- [popgen-vcf-analysis](https://github.com/syrymovmukhammed/popgen-vcf-analysis) — multi-sample VCF → SNP distances, PCA, clustering, NJ tree
