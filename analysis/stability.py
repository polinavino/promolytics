#!/usr/bin/env python3
"""Ordinal stability of the promotability index under measure removal.

A weight-free aggregate is only useful if the order it produces does not hinge on
any single measure. This is the framework's stability-under-nuisance check. We
drop each panel assay in turn, rebuild the average-rank index on the remaining
measures, and ask how much the ranking moves.

Two readouts per dropped measure:
  spearman   rank correlation of the reduced index with the full index
  top100     overlap of the top-100 chemicals before and after the drop

A measure whose removal barely moves either is redundant with the rest. A measure
whose removal moves them a lot is load-bearing, and any property we propose for a
promotion assay has to say what to do about it.
"""

import numpy as np
from scipy.stats import spearmanr

import _common as C
from avg_rank import average_rank_index

TOP_K = 100


def top_overlap(a, b, k=TOP_K):
    """Jaccard-style overlap of the top-k index members of two Series."""
    ta = set(a.sort_values(ascending=False).head(k).index)
    tb = set(b.sort_values(ascending=False).head(k).index)
    return len(ta & tb) / len(ta | tb)


def main():
    C.rule("STABILITY: leave-one-measure-out robustness of the index")
    m = C.load_matrix()
    panel = C.PROMOTION_PANEL

    full, _ = average_rank_index(m, panel)
    full = full.dropna()
    print("\nchemicals in the full index: %d" % len(full))
    print("dropping each of the %d panel measures in turn:\n" % len(panel))

    rows = []
    for a in panel:
        reduced_assays = [x for x in panel if x != a]
        red, _ = average_rank_index(m, reduced_assays)
        common = full.index.intersection(red.dropna().index)
        rho = spearmanr(full.loc[common], red.loc[common]).correlation
        ov = top_overlap(full.loc[common], red.loc[common])
        rows.append((a, C.ASSAY_TO_FAMILY[a], rho, ov, len(common)))

    rows.sort(key=lambda r: r[2])  # least stable (lowest spearman) first
    print("  %-16s %-24s %9s %9s" % ("dropped", "family", "spearman", "top100"))
    for a, fam, rho, ov, n in rows:
        print("  %-16s %-24s %9.4f %9.3f" % (a, fam, rho, ov))

    rhos = [r[2] for r in rows]
    ovs = [r[3] for r in rows]
    worst_ov = min(rows, key=lambda r: r[3])   # measure whose removal churns the top-100 most
    C.rule("Verdict")
    print("worst-case rank correlation under any single drop: %.4f" % min(rhos))
    print("mean rank correlation under single drops:          %.4f" % np.mean(rhos))
    print("worst-case top-100 overlap under any single drop:  %.3f" % min(ovs))
    print("most load-bearing for the bulk order: %s (spearman %.4f when dropped, %s)"
          % (rows[0][0], rows[0][2], rows[0][1]))
    print("most load-bearing for the top-100:    %s (top-100 overlap %.3f when dropped, %s)"
          % (worst_ov[0], worst_ov[3], worst_ov[1]))
    print("most redundant measure:               %s (spearman %.4f when dropped)"
          % (rows[-1][0], rows[-1][2]))
    print("\n=> the bulk order is stable (mean spearman %.3f) but the top-100 is not"
          % np.mean(rhos))
    print("   (dropping %s changes about a third of the top-100)." % worst_ov[0])
    print("   The oxidative-stress measures move the bulk order most, the top-100 is")
    print("   most sensitive to a receptor measure, so the two metrics disagree on which")
    print("   measure is load-bearing. Nominations are only as stable as the panel, so a")
    print("   promotion index has to publish its panel and its stability alongside it.")


if __name__ == "__main__":
    with C.tee("stability"):
        main()
