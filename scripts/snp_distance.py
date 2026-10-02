#!/usr/bin/env python3
"""
snp_distance.py — pairwise SNP distance matrix from the core-SNP VCF.

Writes the square matrix (snp-dists format), a long table, a clustered
heatmap, and — when a truth table exists (sim mode) — the correlation
between measured and true distances, which validates the whole
genome -> alignment -> VCF chain.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from vcf_matrix import read_genotypes, snp_distance


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--truth", help="truth_snp_distance.tsv (sim mode)")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    G, _ = read_genotypes(a.vcf)
    D = snp_distance(G)
    D.to_csv(out / "snp_distance_matrix.tsv", sep="\t")
    iu = np.triu_indices(len(D), 1)
    long = pd.DataFrame({"strain1": D.index[iu[0]], "strain2": D.columns[iu[1]], "snp_distance": D.to_numpy()[iu]})
    long.to_csv(out / "snp_distance_pairs.tsv", sep="\t", index=False)

    order = leaves_list(linkage(squareform(D.to_numpy()), method="average"))
    Do = D.iloc[order, order]
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(Do.to_numpy(), cmap="viridis")
    ax.set_xticks(range(len(Do)), Do.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(Do)), Do.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="SNP distance (core sites)")
    ax.set_title(f"Pairwise SNP distance — {G.shape[0]} strains, {G.shape[1]} core SNPs")
    fig.tight_layout(); fig.savefig(out / "snp_distance_heatmap.png", dpi=130); plt.close(fig)

    stats = {"strains": int(G.shape[0]), "core_snps": int(G.shape[1]),
             "median_pairwise_distance": float(np.median(long["snp_distance"])),
             "max_pairwise_distance": int(long["snp_distance"].max())}
    if a.truth and Path(a.truth).exists():
        t = pd.read_csv(a.truth, sep="\t")
        key = lambda x, y: tuple(sorted((x, y)))
        tmap = {key(x, y): d for x, y, d in t.itertuples(index=False)}
        long["truth"] = [tmap[key(x, y)] for x, y in zip(long.strain1, long.strain2)]
        r = float(np.corrcoef(long["snp_distance"], long["truth"])[0, 1])
        stats.update(pearson_r_vs_truth=round(r, 5),
                     mean_abs_error_vs_truth=round(float((long["snp_distance"] - long["truth"]).abs().mean()), 3))
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(long["truth"], long["snp_distance"], s=10, alpha=.6)
        lim = [0, max(long["truth"].max(), long["snp_distance"].max()) * 1.05]
        ax.plot(lim, lim, ls=":", c="grey")
        ax.set(xlabel="true SNP distance (simulation)", ylabel="measured SNP distance (pipeline)",
               title=f"Pipeline vs truth, r = {r:.4f}")
        fig.tight_layout(); fig.savefig(out / "distance_vs_truth.png", dpi=130); plt.close(fig)
        long.to_csv(out / "snp_distance_pairs.tsv", sep="\t", index=False)
    (out / "distance_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
