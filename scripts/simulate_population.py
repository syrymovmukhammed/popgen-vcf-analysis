#!/usr/bin/env python3
"""
simulate_population.py — evolve a set of bacterial genomes along a known tree.

Tree (default 4 clades x 6 strains), all mutations relative to the reference:

    reference ── clade stem (hundreds of SNPs, different length per clade)
                   ├── subclade a (tens of SNPs) ── strains (private SNPs)
                   └── subclade b (tens of SNPs) ── strains (private SNPs)

Realism added on purpose, because each one is a known source of trouble in
real bacterial population genomics:
  * small indels on every branch          (aligner must handle gaps)
  * accessory-gene loss: 0-2 deletions of 1-4 kb per strain
                                          (missing data -> core-genome filter)
  * one recombination event: a strain imports a 6 kb segment from another
    clade                                 (homoplasy; distorts trees)
Clade names mimic E. coli phylogroups (A, B1, B2, D) — the labels are the
ground truth used to score the clustering.

Outputs (in --outdir):
  genomes/<strain>.fa       one FASTA per strain
  strain_labels.tsv         strain -> clade (truth)
  truth_snp_distance.tsv    exact pairwise SNP distances over shared sequence
  truth_events.tsv          deletions / recombination per strain
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from wgs_sim import read_fasta, write_fasta

CLADES = ["A", "B1", "B2", "D", "E", "F"]
# stem lengths in SNPs per 100 kb: A is closest to the K-12 reference
STEM_PER_100KB = {"A": 60, "B1": 450, "B2": 700, "D": 600, "E": 650, "F": 550}


def mutate(alleles: dict[int, str], ref: str, n: int, rng: np.random.Generator,
           forbid: set[int] | None = None) -> dict[int, str]:
    """Add n new SNPs to an allele map {pos0: base}; returns a new map."""
    out = dict(alleles)
    L = len(ref)
    while n > 0:
        p = int(rng.integers(200, L - 200))
        if forbid and p in forbid:
            continue
        cur = out.get(p, ref[p])
        if cur == "N":
            continue
        out[p] = rng.choice([b for b in "ACGT" if b != cur])
        if out[p] == ref[p]:
            del out[p]
        n -= 1
    return out


def build_genome(ref: str, snps: dict[int, str], indels: list[tuple[int, str, int]],
                 deletions: list[tuple[int, int]]) -> str:
    """Apply SNPs, then small indels and large deletions (all in ref coords)."""
    seq = list(ref)
    for p, b in snps.items():
        seq[p] = b
    edits = [(p, ins, dl) for p, ins, dl in indels] + [(s, "", e - s) for s, e in deletions]
    out, last = [], 0
    for p, ins, dl in sorted(edits):
        if p < last:
            continue                      # overlapping edit -> skip
        out.append("".join(seq[last:p])); out.append(ins); last = p + dl
    out.append("".join(seq[last:]))
    return "".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-clades", type=int, default=4)
    ap.add_argument("--strains-per-clade", type=int, default=6)
    ap.add_argument("--seed", type=int, default=2024)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    _, ref = read_fasta(a.reference)
    L = len(ref)
    scale = L / 100_000
    out = Path(a.outdir); (out / "genomes").mkdir(parents=True, exist_ok=True)

    strains: dict[str, dict] = {}
    for c in CLADES[: a.n_clades]:
        stem = mutate({}, ref, int(STEM_PER_100KB[c] * scale), rng)
        n_a = (a.strains_per_clade + 1) // 2
        for sub, n_in_sub in (("a", n_a), ("b", a.strains_per_clade - n_a)):
            subclade = mutate(stem, ref, int(rng.integers(25, 45) * scale), rng)
            for i in range(n_in_sub):
                name = f"{c}_{sub}{i + 1}"
                snps = mutate(subclade, ref, int(rng.integers(5, 30) * scale), rng)
                indels = []
                for _ in range(int(3 * scale)):
                    p = int(rng.integers(500, L - 500))
                    if rng.random() < 0.5:
                        indels.append((p, "".join(rng.choice(list("ACGT"), int(rng.integers(1, 6)))), 0))
                    else:
                        indels.append((p, "", int(rng.integers(1, 6))))
                dels = []
                for _ in range(int(rng.integers(0, 3))):
                    s = int(rng.integers(2000, L - 6000))
                    dels.append((s, s + int(rng.integers(1000, 4000))))
                strains[name] = {"clade": c, "snps": snps, "indels": indels, "dels": dels, "events": []}

    # one recombination: first strain of clade 2 imports 6 kb from clade 3's stem lineage
    names = sorted(strains)
    if a.n_clades >= 3:
        recip = next(n for n in names if strains[n]["clade"] == CLADES[1])
        donor = next(n for n in names if strains[n]["clade"] == CLADES[2])
        s = int(rng.integers(5000, L - 11000)); e = s + 6000
        r = strains[recip]["snps"]
        for p in [p for p in r if s <= p < e]:
            del r[p]
        for p, b in strains[donor]["snps"].items():
            if s <= p < e:
                r[p] = b
        strains[recip]["events"].append(f"recombination from {donor} at {s + 1}-{e}")

    for n, d in strains.items():
        write_fasta(out / "genomes" / f"{n}.fa", [(n, build_genome(ref, d["snps"], d["indels"], d["dels"]))])
        for s, e in d["dels"]:
            d["events"].append(f"deletion {s + 1}-{e}")

    with open(out / "strain_labels.tsv", "w") as fh:
        fh.write("strain\tgroup\n")
        for n in names:
            fh.write(f"{n}\t{strains[n]['clade']}\n")
    with open(out / "truth_events.tsv", "w") as fh:
        fh.write("strain\tevent\n")
        for n in names:
            for ev in strains[n]["events"]:
                fh.write(f"{n}\t{ev}\n")

    # exact SNP distance over sequence present in both strains (deletions excluded)
    present = {}
    for n in names:
        m = np.ones(L, dtype=bool)
        for s, e in strains[n]["dels"]:
            m[s:e] = False
        present[n] = m
    with open(out / "truth_snp_distance.tsv", "w") as fh:
        fh.write("strain1\tstrain2\tsnp_distance\n")
        for x, y in combinations(names, 2):
            sx, sy = strains[x]["snps"], strains[y]["snps"]
            pos = {p for p in set(sx) | set(sy) if present[x][p] and present[y][p]}
            d = sum(sx.get(p, ref[p]) != sy.get(p, ref[p]) for p in pos)
            fh.write(f"{x}\t{y}\t{d}\n")
    print(f"simulated {len(names)} genomes in {a.n_clades} clades -> {out / 'genomes'}")


if __name__ == "__main__":
    main()
