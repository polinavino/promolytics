#!/usr/bin/env python3
"""Orient the promotability index against an external carcinogenicity anchor.

The average-rank index is sign-free on its own. To read it as promotional
potential we check it against an outcome the concept should track. The anchor is
carcinogenicity from the Lagunin set, and its non-genotoxic subset (carcinogenic
and not mutagenic in Ames), which is the closest public proxy for tumour
promotion.

Two questions:

1. Orientation. Do higher-index chemicals carry carcinogenicity more often than
   lower-index ones? Reported as the rank AUC of the index against the label,
   with a Mann-Whitney p-value.

2. Extremes versus middle. The framework claims a consensus of disagreeing
   measures is trustworthy at its extremes and undecided in its middle. We cut
   the full index distribution into bands and read the label rate in each band.

The anchor is small and noisy, and the script prints every N so the reader can
see how thin it is. Treat the direction as suggestive, not as an effect size.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import _common as C


def auc_rank(scores, labels):
    """AUC via the Mann-Whitney U statistic. Returns (auc, p, n_pos, n_neg)."""
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), len(pos), len(neg)
    u, p = mannwhitneyu(pos, neg, alternative="two-sided")
    auc = u / (len(pos) * len(neg))
    return auc, p, len(pos), len(neg)


def hitcount_bands(df, label_col, buckets):
    """Label rate by number of promotion-pathway hits.

    Quantile bands on the index are degenerate here because most chemicals hit
    zero pathways and pile up at the bottom. The number of promotion hits is the
    interpretable axis: zero hits is the bottom extreme, many hits the top, one
    or two hits the ambiguous middle.
    """
    rows = []
    for name, lo, hi in buckets:
        sel = (df["n_hits"] >= lo) & (df["n_hits"] <= hi)
        sub = df[sel]
        lab = sub[label_col].dropna()
        rate = lab.mean() if len(lab) else float("nan")
        rows.append((name, len(sub), int(lab.notna().sum()),
                     int((lab == 1).sum()), rate))
    return rows


def report_anchor(scored, label_col, title):
    C.rule(title)
    lab = scored.dropna(subset=[label_col])
    auc, p, npos, nneg = auc_rank(lab["promotability"].values, lab[label_col].values)
    print("labeled chemicals with an index: %d  (%d positive / %d negative)"
          % (len(lab), npos, nneg))
    if np.isnan(auc):
        print("not enough labeled chemicals in one class to orient. skipping.")
        return
    print("rank AUC of index vs %s: %.3f  (Mann-Whitney p=%.3g)" % (label_col, auc, p))
    direction = "positively" if auc > 0.5 else "negatively"
    print("=> higher promotability index tracks %s with %s" % (label_col, direction))

    print("\nlabel rate by number of promotion-pathway hits:")
    buckets = [("0 hits", 0, 0), ("1 hit", 1, 1), ("2 hits", 2, 2),
               ("3 hits", 3, 3), ("4+ hits", 4, 99)]
    rows = hitcount_bands(scored, label_col, buckets)
    print("  %-9s %8s %8s %8s %8s" % ("hits", "n", "labeled", "pos", "rate"))
    for name, n, nl, npos_b, rate in rows:
        rate_s = "  n/a" if np.isnan(rate) else "%5.1f%%" % (100 * rate)
        print("  %-9s %8d %8d %8d %8s" % (name, n, nl, npos_b, rate_s))
    print("(zero hits is the bottom extreme, 4+ the top, one or two the middle)")


def main():
    C.rule("ANCHOR: orient promotability against carcinogenicity")
    m = C.load_matrix()
    idx = pd.read_csv(f"{C.DATA}/promotability_index.csv")
    panel = C.PROMOTION_PANEL
    n_hits = (m[panel] == 1).sum(axis=1)
    labels = m[["ikey", "ames", "carcinogen", "ngc"]].assign(n_hits=n_hits.values)
    scored = idx.merge(labels, on="ikey", how="left")

    print("\nchemicals with a promotability index: %d" % len(scored))
    print("  with carcinogen label:  %d" % scored["carcinogen"].notna().sum())
    print("  with ngc label:         %d" % scored["ngc"].notna().sum())

    report_anchor(scored, "carcinogen", "Anchor 1: all carcinogens")
    report_anchor(scored, "ngc", "Anchor 2: non-genotoxic carcinogens (promotion proxy)")

    C.rule("Reading")
    print("The anchor is sparse (tens of positives), so the bands are illustrative.")
    print("What the framework asks is a direction check plus a shape check: the index")
    print("should separate outcomes at its extremes and blur in its middle. Scale-up")
    print("path is ToxCast/invitroDB potencies plus CPDB/IARC carcinogen calls.")


if __name__ == "__main__":
    with C.tee("anchor"):
        main()
