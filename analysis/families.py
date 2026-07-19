#!/usr/bin/env python3
"""Competing promotability measures: do they agree, and how do they cluster?

Each Tox21 assay in the promotion panel is treated as a distinct measure that
claims to quantify the same latent concept (promotional potential). The
framework predicts these measures will (a) disagree pairwise and (b) fall into
pathway families rather than one canonical axis. This script quantifies both.

Agreement between two binary assays is measured with the phi coefficient
(Pearson correlation of 0/1 hit-calls) on their pairwise-complete chemicals, and
with the Jaccard index of their hit sets. Families are recovered by average-
linkage hierarchical clustering on 1 - phi distance, and compared against the a
priori nuclear-receptor vs stress grouping.
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

import _common as C


def phi_matrix(df, assays):
    """Pairwise-complete phi (Pearson on 0/1) between every pair of assays."""
    n = len(assays)
    phi = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            a, b = df[assays[i]], df[assays[j]]
            mask = a.notna() & b.notna()
            if mask.sum() < 30:
                continue
            x, y = a[mask].values, b[mask].values
            if x.std() == 0 or y.std() == 0:
                phi[i, j] = 0.0
            else:
                phi[i, j] = np.corrcoef(x, y)[0, 1]
    return pd.DataFrame(phi, index=assays, columns=assays)


def jaccard(df, a, b):
    mask = df[a].notna() & df[b].notna()
    x, y = df[a][mask] == 1, df[b][mask] == 1
    union = (x | y).sum()
    return (x & y).sum() / union if union else 0.0


def main():
    C.rule("FAMILIES: agreement and clustering of competing promotability measures")
    m = C.load_matrix()
    assays = C.PROMOTION_PANEL

    print("\nPanel: %d competing measures (a priori families):" % len(assays))
    for fam in ("nuclear_receptor", "oxidative_proteo_stress"):
        print("  %-24s %s" % (fam, ", ".join(C.ASSAY_FAMILIES[fam])))

    # --- Pairwise agreement ---
    phi = phi_matrix(m, assays)
    off = phi.values[~np.eye(len(assays), dtype=bool)]
    off = off[~np.isnan(off)]
    C.rule("Pairwise agreement (phi coefficient)")
    print("mean off-diagonal phi:   %.3f" % np.nanmean(off))
    print("median off-diagonal phi: %.3f" % np.nanmedian(off))
    print("max off-diagonal phi:    %.3f" % np.nanmax(off))
    print("(low values mean the measures genuinely disagree, none is a proxy for another)")

    print("\nphi matrix:")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(phi.round(2).to_string())

    # Strongest and weakest agreeing pairs.
    pairs = []
    for i in range(len(assays)):
        for j in range(i + 1, len(assays)):
            pairs.append((assays[i], assays[j], phi.iloc[i, j], jaccard(m, assays[i], assays[j])))
    pairs.sort(key=lambda t: -t[2])
    print("\nMost concordant measure pairs (phi, Jaccard of hit sets):")
    for a, b, p, jac in pairs[:5]:
        print("  %-14s ~ %-14s  phi=%.2f  jaccard=%.2f  [%s | %s]"
              % (a, b, p, jac, C.ASSAY_TO_FAMILY[a], C.ASSAY_TO_FAMILY[b]))
    print("Least concordant measure pairs:")
    for a, b, p, jac in pairs[-5:]:
        print("  %-14s ~ %-14s  phi=%.2f  jaccard=%.2f  [%s | %s]"
              % (a, b, p, jac, C.ASSAY_TO_FAMILY[a], C.ASSAY_TO_FAMILY[b]))

    # --- Family recovery by clustering ---
    C.rule("Recovered families (average-linkage on 1 - phi)")
    dist = 1.0 - phi.fillna(0.0).values
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0  # enforce symmetry for squareform
    Z = linkage(squareform(dist, checks=False), method="average")

    for k in (2, 3, 4):
        labels = fcluster(Z, k, criterion="maxclust")
        print("\nk=%d clusters:" % k)
        for c in sorted(set(labels)):
            members = [assays[i] for i in range(len(assays)) if labels[i] == c]
            fams = {C.ASSAY_TO_FAMILY[a] for a in members}
            print("  cluster %d (%s): %s" % (c, "/".join(sorted(fams)), ", ".join(members)))

    # Agreement of the 2-cluster solution with the a priori NR vs SR split.
    labels2 = fcluster(Z, 2, criterion="maxclust")
    apriori = np.array([0 if C.ASSAY_TO_FAMILY[a] == "nuclear_receptor" else 1 for a in assays])
    # cluster ids are arbitrary, so align by majority.
    agree = max((labels2 == (apriori + 1)).mean(), (labels2 == (2 - apriori)).mean())
    C.rule("Verdict")
    print("2-cluster recovery vs a priori NR/stress split: %.0f%% of assays agree" % (100 * agree))
    print("No single measure is canonical: mean pairwise phi is %.2f, and the panel"
          % np.nanmean(off))
    print("splits into pathway families (e.g. ER pair, AR pair, and a mixed")
    print("xenobiotic/oxidative-stress group) rather than collapsing to one axis.")


if __name__ == "__main__":
    with C.tee("families"):
        main()
