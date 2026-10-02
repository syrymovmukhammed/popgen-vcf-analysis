#!/usr/bin/env python3
"""Collect the key numbers of a run into results/summary.md (also shown as
the GitHub Actions job summary)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    proj = Path(a.project)
    res = proj / "results"
    src = (proj / "data" / "reference.source").read_text().strip()
    d = json.loads((res / "04_distance" / "distance_stats.json").read_text())
    c = json.loads((res / "05_pca_clustering" / "clustering_stats.json").read_text())

    md = ["# Population genomics — run summary", "", f"**Reference:** `{src}`"]
    if src.startswith("SYNTHETIC"):
        md.append("> ⚠️ NCBI was unreachable; this run used the synthetic offline fallback reference.")
    md += ["", "## Data", "",
           f"- strains: **{d['strains']}**, core biallelic SNPs: **{d['core_snps']}**",
           f"- median / max pairwise SNP distance: {d['median_pairwise_distance']:.0f} / {d['max_pairwise_distance']}"]
    if "pearson_r_vs_truth" in d:
        md.append(f"- measured vs true SNP distance: Pearson r = **{d['pearson_r_vs_truth']}**, "
                  f"mean absolute error = {d['mean_abs_error_vs_truth']} SNPs")
    md += ["", "## Structure", "",
           f"- PCA explained variance (PC1–PC5): {', '.join(f'{v:.1%}' for v in c['pca_explained_variance'])}",
           f"- clusters chosen by silhouette: **k = {c['best_k_by_silhouette']}** (silhouette {c['silhouette']})"]
    if "ARI_hierarchical" in c:
        md += [f"- known groups: {c['n_groups_expected']}; agreement with clusters — "
               f"ARI (UPGMA) **{c['ARI_hierarchical']}**, ARI (k-means on PCs) **{c['ARI_kmeans']}**, "
               f"NMI {c['NMI_hierarchical']}", ""]
        ct = pd.read_csv(res / "05_pca_clustering" / "group_vs_cluster.tsv", sep="\t", index_col=0)
        md += ["| group \\ cluster | " + " | ".join(map(str, ct.columns)) + " |",
               "|---|" + "---|" * len(ct.columns)]
        md += [f"| {g} | " + " | ".join(map(str, r)) + " |" for g, r in ct.iterrows()]
    ev = proj / "data" / "truth_events.tsv"
    if ev.exists():
        rec = [l for l in ev.read_text().splitlines() if "recombination" in l]
        if rec:
            md += ["", f"Simulated recombination (expect this strain to sit slightly off its clade): `{rec[0]}`"]
    Path(a.out).write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
