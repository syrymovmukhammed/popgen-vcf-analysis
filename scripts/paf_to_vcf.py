#!/usr/bin/env python3
"""
paf_to_vcf.py — whole-genome alignments -> multi-sample haploid SNP VCF.

Input: one PAF per strain from `minimap2 -c --cs -x asm20 reference.fa strain.fa`.
For every strain the `cs` difference string of each primary alignment block
is walked along the reference:

    :120      120 identical bases          -> covered, reference allele
    *ag       reference a, strain g        -> covered, SNP
    -ttg      bases missing in the strain  -> not covered (missing)
    +ca       bases inserted in the strain -> ignored (SNP-only analysis)

A reference position covered by >1 block of the same strain (repeat / paralog)
is set to missing, so ambiguous sites never become fake SNPs. A site is
emitted when at least one strain carries a non-reference base; each strain's
genotype is ALT index, 0 (reference) or '.' (not covered). This is the logic
behind snippy-core / ParSNP-style core-genome SNP tables, in ~100 lines.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

from wgs_sim import read_fasta

CS_TOKEN = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][a-z]+|~[a-z]{2}\d+[a-z]{2})")
B2I = {b: i for i, b in enumerate("ACGT")}


def walk_paf(path: Path, ref_name: str, L: int, min_len: int):
    """Return (coverage count per ref position, {pos0: alt_base})."""
    cov = np.zeros(L, dtype=np.int16)
    snps: dict[int, str] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[5] != ref_name or int(f[10]) < min_len:
                continue
            tags = dict(t.split(":", 2)[::2] for t in f[12:])
            if tags.get("tp", "P") != "P" or "cs" not in tags:
                continue
            t = int(f[7])
            for tok in CS_TOKEN.findall(tags["cs"]):
                op = tok[0]
                if op == ":":
                    n = int(tok[1:]); cov[t:t + n] += 1; t += n
                elif op == "*":
                    cov[t] += 1; snps[t] = tok[2].upper(); t += 1
                elif op == "-":
                    t += len(tok) - 1
                # '+' (insertion in strain) does not advance the reference
    return cov, snps


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--paf", nargs="+", required=True, help="<strain>.paf files")
    ap.add_argument("--min-aln-len", type=int, default=1000)
    ap.add_argument("--out", required=True, help="uncompressed VCF")
    a = ap.parse_args()

    ref_name, ref = read_fasta(a.reference)
    L = len(ref)
    names, calls = [], []
    for p in sorted(a.paf):
        name = Path(p).stem
        cov, snps = walk_paf(Path(p), ref_name, L, a.min_aln_len)
        ok = cov == 1
        names.append(name)
        calls.append((ok, {k: v for k, v in snps.items() if ok[k]}))
        print(f"{name}: {ok.mean():.1%} of reference aligned uniquely, {len(calls[-1][1])} SNPs",
              file=sys.stderr)

    sites = sorted({p for _, s in calls for p in s})
    with open(a.out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=paf_to_vcf.py (minimap2 cs-tag SNPs)\n")
        fh.write(f"##contig=<ID={ref_name},length={L}>\n")
        fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Alt allele count">\n')
        fh.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of called strains">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for pos in sites:
            r = ref[pos]
            gts, alts = [], []
            for ok, s in calls:
                if not ok[pos]:
                    gts.append(".")
                elif pos in s and s[pos] != r:
                    if s[pos] not in alts:
                        alts.append(s[pos])
                    gts.append(str(alts.index(s[pos]) + 1))
                else:
                    gts.append("0")
            if not alts:
                continue
            ac = [gts.count(str(i + 1)) for i in range(len(alts))]
            an = sum(g != "." for g in gts)
            fh.write(f"{ref_name}\t{pos + 1}\t.\t{r}\t{','.join(alts)}\t.\tPASS\t"
                     f"AC={','.join(map(str, ac))};AN={an}\tGT\t" + "\t".join(gts) + "\n")
    print(f"{len(sites)} variable sites across {len(names)} strains -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
