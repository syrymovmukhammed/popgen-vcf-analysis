#!/usr/bin/env python3
"""
pca_clustering.py — population structure from the core-SNP VCF.

1. PCA of the strain x SNP genotype matrix (missing -> site mean).
2. Hierarchical clustering (UPGMA / average linkage) of the SNP-distance
   matrix; the number of clusters k is chosen by the silhouette score, not
   fixed by hand. K-means on the leading PCs is run as a second, independent
   method with the same k.
3. Neighbour-joining tree (Biopython) from the same distance matrix,
   midpoint-rooted, written as Newick and drawn as PNG.
4. Agreement of both clusterings with the known groups (simulated clades or
   literature phylogroups): Adjusted Rand Index and NMI.
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
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

from vcf_matrix import read_genotypes, snp_distance

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#7f7f7f"]


def nj_tree(D: pd.DataFrame):
    names = list(D.index)
    lower = [[float(D.iloc[i, j]) for j in range(i + 1)] for i in range(len(names))]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, lower))
    for c in tree.get_nonterminals():
        c.name = None             # drop "Inner1..." labels for a clean Newick
    tree.root_at_midpoint()
    tree.ladderize()
    return tree


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--labels", help="strain<TAB>group table (truth / expected)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-k", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2024)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    G, _ = read_genotypes(a.vcf)
    strains = list(G.index)
    X = G.to_numpy()
    X = np.where(np.isnan(X), np.nanmean(X, axis=0), X)
    D = snp_distance(G)
    groups = None
    if a.labels and Path(a.labels).exists():
        lab = pd.read_csv(a.labels, sep="\t").set_index("strain")["group"]
        groups = lab.reindex(strains).fillna("NA").to_numpy()

    # ---- PCA
    pca = PCA(n_components=min(10, len(strains) - 1), random_state=a.seed)
    pcs = pca.fit_transform(X)
    ev = pca.explained_variance_ratio_

    # ---- hierarchical clustering, k by silhouette
    Z = linkage(squareform(D.to_numpy()), method="average")
    sil = {}
    for k in range(2, min(a.max_k, len(strains) - 1) + 1):
        lab_k = fcluster(Z, k, criterion="maxclust")
        if len(set(lab_k)) > 1:
            sil[k] = silhouette_score(D.to_numpy(), lab_k, metric="precomputed")
    best_k = max(sil, key=sil.get)
    hc = fcluster(Z, best_k, criterion="maxclust")
    km = KMeans(n_clusters=best_k, n_init=20, random_state=a.seed).fit_predict(pcs[:, :min(5, pcs.shape[1])])

    tbl = pd.DataFrame({"strain": strains, "group": groups if groups is not None else "NA",
                        "hc_cluster": hc, "kmeans_cluster": km + 1,
                        "PC1": pcs[:, 0], "PC2": pcs[:, 1]})
    tbl.to_csv(out / "clusters.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame({"k": list(sil), "silhouette": list(sil.values())}).to_csv(
        out / "silhouette_by_k.tsv", sep="\t", index=False, float_format="%.4f")

    stats = {"strains": len(strains), "core_snps": int(G.shape[1]),
             "pca_explained_variance": [round(float(v), 4) for v in ev[:5]],
             "best_k_by_silhouette": int(best_k), "silhouette": round(float(sil[best_k]), 4)}
    if groups is not None:
        stats.update(
            n_groups_expected=int(len(set(groups))),
            ARI_hierarchical=round(adjusted_rand_score(groups, hc), 4),
            ARI_kmeans=round(adjusted_rand_score(groups, km), 4),
            NMI_hierarchical=round(normalized_mutual_info_score(groups, hc), 4))
        ct = pd.crosstab(pd.Series(groups, name="group"), pd.Series(hc, name="hc_cluster"))
        ct.to_csv(out / "group_vs_cluster.tsv", sep="\t")

    # ---- NJ tree
    tree = nj_tree(D)
    Phylo.write(tree, out / "nj_tree.nwk", "newick")
    gcol = {}
    if groups is not None:
        for i, g in enumerate(sorted(set(groups))):
            gcol[g] = PALETTE[i % len(PALETTE)]
    fig = plt.figure(figsize=(7, max(5, 0.28 * len(strains))))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False,
               label_colors={s: gcol.get(g, "k") for s, g in zip(strains, tbl["group"])})
    ax.set_title("Neighbour-joining tree (SNP distance, midpoint-rooted)")
    fig.tight_layout(); fig.savefig(out / "nj_tree.png", dpi=130); plt.close(fig)

    # ---- PCA plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for c in sorted(set(hc)):
        m = hc == c
        ax[0].scatter(pcs[m, 0], pcs[m, 1], s=55, color=PALETTE[(c - 1) % len(PALETTE)],
                      label=f"cluster {c}", edgecolor="k", lw=.4)
    if groups is not None:            # one label per known group, at its centroid
        for g in sorted(set(groups)):
            m = groups == g
            ax[0].annotate(g, (pcs[m, 0].mean(), pcs[m, 1].mean()), fontsize=10, weight="bold",
                           xytext=(8, 8), textcoords="offset points")
    else:
        for i, s in enumerate(strains):
            ax[0].annotate(s, (pcs[i, 0], pcs[i, 1]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[0].set(xlabel=f"PC1 ({ev[0]:.1%})", ylabel=f"PC2 ({ev[1]:.1%})",
              title=f"PCA of core SNPs — colour = cluster (k={best_k}), text = known group")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    dendrogram(Z, labels=strains, ax=ax[1], leaf_font_size=7, color_threshold=Z[-best_k + 1, 2])
    ax[1].set(title="UPGMA dendrogram (SNP distance)", ylabel="SNP distance")
    fig.tight_layout(); fig.savefig(out / "pca_clusters.png", dpi=130); plt.close(fig)

    (out / "clustering_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
