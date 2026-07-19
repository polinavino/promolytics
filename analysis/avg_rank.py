#!/usr/bin/env python3
"""Consensus poset and the average-rank promotability index.

The framework's canonical aggregate of many disagreeing measures is weight-free.
Two objects are built here.

1. Consensus poset. Chemical A ranks at least as high as chemical B under the
   consensus only if A hits a superset of the promotion assays that B hits. This
   is Pareto dominance on the binary hit-vectors. Where two chemicals hit
   different, non-nested sets of assays the consensus abstains, so they are left
   incomparable. The poset records exactly where all measures agree on order.

2. Average-rank promotability index. Each measure is turned into a percentile
   rank of chemicals (a hit outranks a non-hit), and each chemical's index is the
   mean of its percentile ranks over the assays on which it was tested. This is a
   complete order, and it is the weight-free aggregate the framework recommends.

The interesting quantity is the gap between the two. The average-rank index
imposes an order on many pairs the poset leaves incomparable. Those forced
orderings are the near-ties, and they mark where the concept, not the data,
is doing the deciding.
"""

from itertools import combinations
from math import comb

import numpy as np
import pandas as pd

import _common as C

MIN_TESTED = 5     # a chemical needs at least this many panel measures for a stable index
NEAR_TIE_EPS = 0.01  # incomparable pairs whose index differs by less than this are near-ties


def percentile_ranks(m, assays):
    """Per-assay percentile rank of chemicals (hit outranks non-hit)."""
    ranks = pd.DataFrame(index=m.index)
    for a in assays:
        col = m[a]
        # rank only chemicals tested on this assay, average ties, scale to (0,1]
        ranks[a] = col.rank(method="average", pct=True)
    return ranks


def average_rank_index(m, assays, min_tested=MIN_TESTED):
    """Mean percentile rank across tested measures, for adequately-tested chemicals."""
    ranks = percentile_ranks(m, assays)
    tested = m[assays].notna().sum(axis=1)
    idx = ranks.mean(axis=1, skipna=True)
    idx[tested < min_tested] = np.nan
    return idx, tested


def poset_stats(patterns_counts):
    """Comparability and maximal elements over distinct hit-patterns.

    patterns_counts: dict frozenset(hit assays) -> number of chemicals.
    Returns a dict of summary statistics computed at the chemical-pair level.
    """
    pats = list(patterns_counts)
    counts = np.array([patterns_counts[p] for p in pats])
    n = int(counts.sum())
    total_pairs = comb(n, 2)

    # tied pairs (same pattern) are comparable ties
    tied = int(sum(comb(c, 2) for c in counts))

    comparable_cross = 0
    for i, j in combinations(range(len(pats)), 2):
        pi, pj = pats[i], pats[j]
        if pi <= pj or pj <= pi:  # subset relation => Pareto comparable
            comparable_cross += int(counts[i] * counts[j])

    comparable = tied + comparable_cross
    # maximal patterns: not a strict subset of any other pattern
    maximal = [p for p in pats if not any(p < q for q in pats)]
    maximal_chems = int(sum(patterns_counts[p] for p in maximal))
    return {
        "n": n,
        "total_pairs": total_pairs,
        "comparable": comparable,
        "tied": tied,
        "incomparable": total_pairs - comparable,
        "frac_comparable": comparable / total_pairs if total_pairs else float("nan"),
        "n_patterns": len(pats),
        "maximal_patterns": maximal,
        "maximal_chems": maximal_chems,
    }


def main():
    C.rule("AVERAGE RANK: consensus poset + weight-free promotability index")
    m = C.load_matrix()
    panel = C.PROMOTION_PANEL

    idx, tested = average_rank_index(m, panel)
    m = m.assign(promotability=idx, n_tested=tested)
    scored = m[m["promotability"].notna()].copy()
    print("\nchemicals with a stable index (>= %d measures tested): %d"
          % (MIN_TESTED, len(scored)))
    print("promotability index range: %.3f .. %.3f, mean %.3f"
          % (scored["promotability"].min(), scored["promotability"].max(),
             scored["promotability"].mean()))

    # --- Consensus poset on complete cases ---
    C.rule("Consensus poset (Pareto dominance on hit-vectors, complete cases)")
    complete = m[m[panel].notna().all(axis=1)].copy()
    hit = complete[panel].astype(int)
    patterns = {}
    for _, row in hit.iterrows():
        key = frozenset(a for a in panel if row[a] == 1)
        patterns[key] = patterns.get(key, 0) + 1

    stats_all = poset_stats(patterns)
    print("complete-case chemicals:        %d" % stats_all["n"])
    print("distinct hit-patterns:          %d" % stats_all["n_patterns"])
    print("fraction of pairs comparable:   %.3f" % stats_all["frac_comparable"])
    print("  (of which tied, same pattern: %.3f)" % (stats_all["tied"] / stats_all["total_pairs"]))
    print("incomparable pairs:             %d of %d"
          % (stats_all["incomparable"], stats_all["total_pairs"]))
    print("maximal (undominated) patterns: %d, covering %d chemicals"
          % (len(stats_all["maximal_patterns"]), stats_all["maximal_chems"]))

    # The empty pattern (no promotion hits) is the unique bottom and inflates
    # comparability. Recompute among chemicals with at least one promotion hit.
    active = {p: c for p, c in patterns.items() if len(p) > 0}
    stats_active = poset_stats(active)
    print("\nAmong active chemicals only (>= 1 promotion hit):")
    print("  active chemicals:             %d" % stats_active["n"])
    print("  fraction of pairs comparable: %.3f" % stats_active["frac_comparable"])
    print("  incomparable pairs:           %d of %d"
          % (stats_active["incomparable"], stats_active["total_pairs"]))
    print("  => the consensus abstains on %.0f%% of active pairs"
          % (100 * (1 - stats_active["frac_comparable"])))

    print("\nMaximal (undominated) promotion profiles, broadest first:")
    for p in sorted(stats_active["maximal_patterns"], key=lambda s: -len(s))[:8]:
        print("  %2d hits: %s  (%d chemicals)"
              % (len(p), ", ".join(sorted(p)), active[p]))

    # --- Does the average-rank order respect the poset? ---
    C.rule("Average-rank order vs the poset")
    # For each comparable (strict) pattern pair, the aggregate should order the
    # dominating pattern's chemicals above the dominated pattern's. Check at the
    # pattern level using the mean index per pattern.
    pat_idx = {}
    comp_scored = complete[complete["promotability"].notna()]
    hit_s = comp_scored[panel].astype(int)
    for i, (_, row) in enumerate(hit_s.iterrows()):
        key = frozenset(a for a in panel if row[a] == 1)
        pat_idx.setdefault(key, []).append(comp_scored["promotability"].iloc[i])
    pat_mean = {p: float(np.mean(v)) for p, v in pat_idx.items()}
    inversions = 0
    checked = 0
    for p, q in combinations(pat_mean, 2):
        if p < q:
            checked += 1
            if pat_mean[p] > pat_mean[q]:
                inversions += 1
        elif q < p:
            checked += 1
            if pat_mean[q] > pat_mean[p]:
                inversions += 1
    print("strictly comparable pattern pairs checked: %d" % checked)
    print("order inversions (aggregate contradicts consensus): %d" % inversions)
    print("=> the average-rank index is a linear extension of the consensus poset"
          if inversions == 0 else
          "=> the average-rank index contradicts the consensus on %d pattern pairs" % inversions)

    # --- Near-ties: incomparable pairs the aggregate still orders ---
    C.rule("Near-ties (consensus abstains, aggregate still decides)")
    act = comp_scored[comp_scored[panel].astype(int).sum(axis=1) > 0].copy()
    act = act.reset_index(drop=True)
    act_hit = [frozenset(a for a in panel if act[a].iloc[i] == 1) for i in range(len(act))]
    act_val = act["promotability"].values
    # sample near-ties for reporting (full O(n^2) over active chemicals is fine here)
    near = []
    for i in range(len(act)):
        for j in range(i + 1, len(act)):
            if act_hit[i] <= act_hit[j] or act_hit[j] <= act_hit[i]:
                continue  # comparable, not a near-tie
            d = abs(act_val[i] - act_val[j])
            if d < NEAR_TIE_EPS:
                near.append((d, i, j))
    near.sort()
    print("active chemicals compared:           %d" % len(act))
    total_active_pairs = comb(len(act), 2)
    print("incomparable pairs within eps=%.2f:  %d" % (NEAR_TIE_EPS, len(near)))
    print("  as a share of all active pairs:    %.3f" % (len(near) / total_active_pairs))
    print("\nexample near-ties (aggregate splits them, consensus does not):")
    for d, i, j in near[:6]:
        print("  d=%.4f  {%s}  vs  {%s}"
              % (d, ", ".join(sorted(act_hit[i])), ", ".join(sorted(act_hit[j]))))

    # --- Top of the index ---
    C.rule("Top promotability (average-rank index)")
    top = scored.sort_values("promotability", ascending=False).head(12)
    for _, r in top.iterrows():
        hits = [a for a in panel if r[a] == 1]
        print("  %.3f  %s  hits=%d {%s}"
              % (r["promotability"], r["ikey"], len(hits), ", ".join(hits)))

    # persist the index for downstream anchor analysis
    out = scored[["ikey", "smiles", "promotability", "n_tested"]]
    out.to_csv(f"{C.DATA}/promotability_index.csv", index=False)
    print("\nwrote data/promotability_index.csv (%d chemicals)" % len(out))


if __name__ == "__main__":
    with C.tee("avg_rank"):
        main()
