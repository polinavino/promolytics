# Properties of a well-formed promotion measure

The field wants a promotion assay, the analogue of the Ames test for the
promotion step. Before benchmarking candidates against each other, it is worth
stating what a candidate has to satisfy to count as well-formed. Each property
below is tied to a concrete failure mode seen in the audit in this repo, so none
of them is abstract. The numbers cited are the current tracked `outputs/`.

These are properties, not a ranking. A measure can be useful without meeting all
of them, but a measure that meets none should not be read as a promotability
score.

## 1. Declared panel and mechanism

A promotion score has to name the measures it aggregates and the pathway each one
reads. Promotional potential is not one axis. Across the 10-measure panel the mean
pairwise agreement (phi) is only 0.19, and the measures cluster into pathway
families (a receptor-ligand group, an androgen group, a mixed
xenobiotic/oxidative-stress group) rather than collapsing onto a single
dimension. A score reported without its panel hides which mechanisms it actually
covers. Mapped onto the key characteristics of carcinogens (ext 3), the panel
reads only 4 of 10 characteristics, with no readout for chronic inflammation,
so the declared panel also has to declare what it does not cover.

*Failure mode it prevents:* a single-assay "promotion score" presented as if it
were the whole concept.

## 2. Reliability and domain gate

A measure has to say for which chemicals it is defined and refuse to score the
rest. In the raw panel, per-assay missingness runs from 7% to 26%, and only 3,013
of 7,492 chemicals are tested on all 10 measures. The index in this repo is
reported only for chemicals tested on at least 5 measures, and this gate should
be explicit in any promotion score rather than silently imputed to zero.

*Failure mode it prevents:* untested chemicals scored as inactive by default,
which understates promotional potential for exactly the compounds no one measured.

## 3. Order faithfulness to the consensus

Where every measure agrees on the order of two chemicals, the aggregate must not
reverse it. The average-rank index here is a linear extension of the consensus
Pareto poset, with zero order inversions against it. Any proposed promotion score
should be checked against the poset the same way, because a score that contradicts
unanimous evidence is not aggregating it.

*Failure mode it prevents:* a weighted score whose weights flip an order that all
constituent measures agreed on.

## 4. Honesty about the undecided middle

A promotion score has to distinguish the region where its measures agree from the
region where they do not. The consensus poset leaves 72% of active-chemical pairs
incomparable, and the average-rank index forces an order on 40% of active pairs
that are genuine near-ties (equal pathway breadth, different pathways). A
well-formed score reports a comparability or confidence flag next to each pair, so
a forced ordering in the middle is not mistaken for a real difference. Grading the
measures by potency (ext 1) is one concrete cure: qHTS AC50 values give a distinct
value to 87% of the chemicals the binary index left tied.

*Failure mode it prevents:* reading a hairline gap in a total order as a
meaningful ranking of two chemicals the evidence cannot separate.

## 5. Ordinal stability under nuisance

The order a score produces should not hinge on any one measure. Dropping any
single panel measure keeps the bulk order almost fixed (mean Spearman 0.96, with
the oxidative-stress measures moving it most), but the top-100 nominations churn by
up to a third, and there a receptor measure is the load-bearing one (dropping
NR-AhR or NR-ER-LBD changes about a third of the top-100). A promotion score used
to nominate chemicals must publish this stability, because the nominations are
only as stable as the panel behind them.

*Failure mode it prevents:* a shortlist that quietly depends on the inclusion of
one favoured assay.

## 6. Monotonicity against an outcome anchor

Higher score should track a promotion-relevant outcome, at least at the extremes.
The current anchor is too thin to confirm this (110 carcinogen-labeled chemicals,
21 positive, rank AUC 0.55 and not significant, and only 8 non-genotoxic
carcinogens), which is itself the point: a promotion score has to declare its
anchor and its power, and a score with no outcome check is unfalsifiable. The
anchor also has to be the right one. Expanding it tenfold with IARC and NTP calls
(ext 2) does not help, because that list is selection-biased (about 57% positive
even at zero promotion hits) and about half of it is genotoxic initiation the panel
does not read, so the index stays flat against it. Restricting to the non-genotoxic
(Ames-negative) subset does not fix the overall separation (AUC 0.51) but the top
of the index is enriched among active chemicals, which is the honest, partial win.

*Failure mode it prevents:* a plausible-looking index tied to an anchor that is
either absent or measuring the wrong thing.

## 7. Honest generalization estimate

When a promotion score is learned from structure, its reported accuracy has to
come from a split that holds out whole scaffolds. Under a random split the same
models score 0.036 to 0.084 AUC higher than under a scaffold-held-out split across
NR-AhR, SR-ARE and carcinogenicity. The honest number is the scaffold one, and a
promotion predictor reported on random folds overstates how it will perform on new
chemistry.

*Failure mode it prevents:* a benchmark leaderboard that rewards memorizing known
scaffolds instead of generalizing to new ones.
