"""Read a haploid multi-sample VCF into a genotype matrix (strains x sites).
Values: 0 = reference, 1 = alternate, NaN = missing. Only biallelic records
are expected (step 03 filters to them)."""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd


def read_genotypes(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (G, sites): G is a strains x sites float DataFrame."""
    opener = gzip.open if path.endswith(".gz") else open
    samples, rows, sites = [], [], []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            gt = [x.split(":")[0] for x in f[9:]]
            rows.append([np.nan if g in (".", "./.") else float(g != "0") for g in gt])
            sites.append((f[0], int(f[1]), f[3], f[4]))
    G = pd.DataFrame(np.array(rows, dtype=float).T if rows else np.empty((len(samples), 0)),
                     index=samples, columns=[f"{c}:{p}" for c, p, _, _ in sites])
    return G, pd.DataFrame(sites, columns=["chrom", "pos", "ref", "alt"])


def snp_distance(G: pd.DataFrame) -> pd.DataFrame:
    """Pairwise count of differing sites, using only sites called in both strains
    (the same definition as `snp-dists`)."""
    X = G.to_numpy()
    called = ~np.isnan(X)
    Xa = np.nan_to_num(X, nan=0.0)
    n = X.shape[0]
    D = np.zeros((n, n), dtype=int)
    for i in range(n):
        both = called[i] & called
        diff = (Xa[i] != Xa) & both
        D[i] = diff.sum(axis=1)
    return pd.DataFrame(D, index=G.index, columns=G.index)
