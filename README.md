# promolytics

A measurement audit of chemical **promotional potential**, built on public
toxicology data. The question is the one the PROMINENT consortium has posed for
the wet lab, asked here on the computational side: several mechanistic assays each
claim to quantify a chemical's ability to drive the promotion step of
carcinogenesis, they disagree, and no single one is canonical. This repo maps
where they agree, builds a weight-free promotability index from them, checks that
index against a carcinogenicity outcome, and shows how honest evaluation changes
the picture. See [HANDOFF.md](HANDOFF.md) for the scoping and rationale.

This is the first project described in the handoff (the promotability measurement
audit, plus the properties note). It stands independent of any de-novo compound
design, which is out of near-term reach.

## In plain language

Cancer usually needs two things. First a cell picks up a mutation (initiation).
Then something pushes that mutated cell to actually grow into a tumour (promotion).
Some chemicals act as promoters. A major open goal, pursued by the PROMINENT
consortium, is a lab test that measures how strongly a chemical promotes cancer, the
way the classic Ames test measures how strongly a chemical mutates DNA. No such
promotion test exists yet.

The problem is that there is no single agreed measurement. Many existing lab assays
each capture one piece of promotion-like activity, such as inflammation, hormone
receptor activation, oxidative stress, or extra cell growth. They do not agree with
each other, and there is no accepted way to fold them into one promotion score.

This repository treats that as a measurement problem and does four things on public
data:
- it measures how much the assays disagree, and finds they mostly do, so none of
  them is "the" promotion measure.
- it combines them into one score without inventing arbitrary weights, and it is
  honest about the many chemicals the assays simply cannot rank against each other.
- it checks whether that score lines up with which chemicals are actually known
  carcinogens, and reports plainly that the current public labels are too thin and
  too biased toward mutation-driven carcinogens to settle the question yet.
- it checks whether a model trained to predict promotion is as good as it looks, and
  shows that its accuracy drops once it is tested on chemistry it has not seen.

Two side results are worth noting. The public assay panel barely measures chronic
inflammation, which is central to promotion, so there is a concrete data gap. And a
search for compounds that reverse a promotion-like gene pattern turns up drug classes
already studied for cancer prevention, while the reverse search correctly flags the
textbook tumour promoters as a sanity check.

The point is not a finished promotion test. It is an honest map of what today's
public data can and cannot support, and of what a good promotion measure would have
to satisfy. That map is meant to complement the consortium's laboratory work.

## What promotional potential is here

Classical carcinogenesis runs initiation (a genotoxic mutation), promotion
(non-genotoxic expansion of initiated cells, driven by receptor signalling,
inflammation, oxidative stress and proliferation), then progression. Promolytics
are the proposed class of agents that block or reverse promotion. There is no
public "promotion assay", so we proxy promotional potential with the mechanistic
readouts that promotion is known to run through, and we proxy the outcome with
non-genotoxic carcinogenicity (carcinogenic but not mutagenic).

## Data provenance

All data is public and needs no account. The core three come from the
Therapeutics Data Commons mirror on the Harvard Dataverse via
[analysis/download_data.py](analysis/download_data.py). The extensions add three
more public sources, fetched and cached by their own scripts.

| dataset | role | source | size |
|---|---|---|---|
| Tox21 (12 assays) | competing promotability measures | TDC / Dataverse | 7,831 chemicals |
| Ames (Hansen) | mutagenicity / genotoxicity, to subtract initiation | TDC / Dataverse | 7,278 chemicals |
| Carcinogens (Lagunin) | carcinogenicity outcome anchor | TDC / Dataverse | 280 chemicals |
| Tox21 qHTS AC50 | potency for the same pathways (ext 1) | PubChem BioAssay | 10 assays queried, 9 usable |
| IARC + NTP carcinogen calls | larger outcome anchor (ext 2) | PubChem annotations | ~1,900 chemicals |
| L1000 connectivity | in-silico promolytic screen (ext 5) | L1000CDS2 (MaayanLab) | API |

Chemicals are harmonized on the InChIKey connectivity block (salt-stripped,
stereochemistry-insensitive), which is the right granularity for joining assay
records to bioassay labels across sources. PubChem records join through CID to the
same InChIKey block.

## The measure-selection rule (stated up front)

The 12 Tox21 endpoints split into three mechanistic families. The **promotability
panel** is the 10 endpoints that read promotion-relevant key characteristics of
carcinogens. The two DNA-damage-response readouts are genotoxicity-leaning
(initiation, not promotion) and are held out of the panel, though kept in the
matrix so the audit can show they behave differently.

- **nuclear receptor** (in panel): NR-AhR, NR-PPAR-gamma, NR-ER, NR-ER-LBD, NR-AR, NR-AR-LBD, NR-Aromatase
- **oxidative / proteotoxic stress** (in panel): SR-ARE, SR-HSE, SR-MMP
- **DNA damage** (held out): SR-p53, SR-ATAD5

## Pipeline

Run everything with [analysis/run_all.sh](analysis/run_all.sh). It finds a local
interpreter that has numpy, pandas, scipy, scikit-learn and rdkit (the miniforge
`molml` env works), or set `PYTHON` yourself. Every script re-tees its console
output into `outputs/<script>.txt`, so the tracked `outputs/` folder is always the
current result of the current code on the current data. Outputs carry no
timestamps and every random step is seeded, so a diff of `outputs/` shows only
real changes in results.

| script | what it does | output |
|---|---|---|
| [download_data.py](analysis/download_data.py) | fetch the three public datasets | `outputs/download_data.txt` |
| [build_matrix.py](analysis/build_matrix.py) | harmonize into a chemical x assay matrix plus labels | `outputs/build_matrix.txt` |
| [families.py](analysis/families.py) | agreement and clustering of the competing measures | `outputs/families.txt` |
| [avg_rank.py](analysis/avg_rank.py) | consensus poset and the average-rank promotability index | `outputs/avg_rank.txt` |
| [stability.py](analysis/stability.py) | leave-one-measure-out robustness of the index | `outputs/stability.txt` |
| [anchor.py](analysis/anchor.py) | orient the index against carcinogenicity | `outputs/anchor.txt` |
| [honest_eval.py](analysis/honest_eval.py) | random vs scaffold-held-out cross-validation | `outputs/honest_eval.txt` |
| [ext_potency.py](analysis/ext_potency.py) | potency-weighted index from Tox21 qHTS AC50 | `outputs/ext_potency.txt` |
| [ext_anchor.py](analysis/ext_anchor.py) | larger IARC/NTP anchor and its non-genotoxic subset | `outputs/ext_anchor.txt` |
| [ext_kcc.py](analysis/ext_kcc.py) | panel coverage of the key characteristics of carcinogens | `outputs/ext_kcc.txt` |
| [ext_lincs.py](analysis/ext_lincs.py) | in-silico promolytic screen by L1000 connectivity | `outputs/ext_lincs.txt` |

The extensions call PubChem and MaayanLab over the network and cache into
`data/raw`, so they are slow on a cold cache and near-instant once cached. Run
`ext_potency.py` before `ext_anchor.py`: they share the CID to InChIKey cache and
must not run at the same time. [analysis/pubchem_util.py](analysis/pubchem_util.py)
holds the shared PubChem retry and mapping helpers.

## Results at a glance

Numbers below are the current tracked outputs. The files in `outputs/` are the
source of truth.

- **The measures genuinely disagree.** Mean pairwise phi across the 10-measure
  panel is 0.19. The strongest agreement is within receptor sub-pairs (NR-AR with
  NR-AR-LBD at 0.56, NR-ER with NR-ER-LBD at 0.53). Clustering recovers pathway
  families rather than one axis, and only 70% of assays fall on the a priori
  nuclear-receptor versus stress split, because AhR, PPAR-gamma and the
  oxidative-stress assays co-activate.
- **The consensus abstains more than it decides.** On complete cases (3,013
  chemicals fully tested), Pareto dominance on the hit-vectors leaves 72% of pairs
  of active chemicals incomparable. The average-rank index is a faithful linear
  extension of the poset (zero order inversions), yet it forces an order on 40% of
  active pairs that are near-ties (same number of pathway hits, different
  pathways). That 40% is where the concept, not the data, is deciding.
- **The index is stable in bulk but not at the top.** Dropping any single measure
  keeps the overall order almost fixed (mean Spearman 0.96), and the oxidative-
  stress measures move it most (dropping SR-ARE gives the lowest correlation,
  0.86). The top-100 nominations are far less stable, and there a receptor measure
  is the load-bearing one: dropping NR-AhR or NR-ER-LBD changes about a third of
  the top-100 (overlap 0.64). The two metrics disagree on which measure matters
  most, so nominations are only as stable as the exact panel that produced them.
- **The index does not yet validate against a carcinogenicity outcome, and that is
  an honest null.** The only positive-looking number, a rank AUC of 0.55 against all
  carcinogens, sits on the least valid anchor, which is contaminated by genotoxic
  initiators the promotion panel is not built to detect, and it is not significant.
  The cleaner non-genotoxic anchors do not orient the index either (rank AUC 0.335
  on the 8 overlapping non-genotoxic carcinogens, and 0.514 on the expanded set in
  ext 2). Only 110 index-scored chemicals carry any carcinogen label (21 positive).
  This is the crude-and-thin-anchor limitation the handoff predicted. It does not
  affect the index or the poset, which stand on their own, but outcome validation
  waits on a cleaner, promotion-specific anchor.
- **Honest evaluation costs accuracy.** Predicting each endpoint from structure
  (ECFP4, random forest) loses AUC under a scaffold-held-out split versus a random
  split: NR-AhR 0.891 to 0.806, SR-ARE 0.793 to 0.719, carcinogen 0.867 to 0.831.
  A promotion model reported on random folds overstates how well it calls unseen
  chemistry. Report the scaffold number.

The properties a well-formed promotion assay should satisfy, each tied to one of
these failure modes, are in [PROPERTIES.md](PROPERTIES.md).

## Extension results

Four of the five extensions in the original plan are built and run on real public
data. Numbers are the current tracked outputs.

- **Potency grades what binary hit-calls tie (ext 1).** Replacing hit-calls with
  Tox21 qHTS AC50 potencies (9 pathways, 8,281 chemicals) gives an index that
  agrees with the binary one in bulk (Spearman 0.86) but assigns a distinct value
  to 87% of the chemicals the binary index left tied. Potency is the natural cure
  for the near-tie problem the binary consensus creates.
- **A bigger anchor is not a better anchor (ext 2).** Pulling IARC and NTP
  carcinogen calls lifts the labelled-and-scored overlap from 110 to 668, but the
  index does not orient against it (rank AUC 0.49, label rate flat across
  promotion-hit counts). The list is selection-biased, since chemicals are
  evaluated because they are suspected carcinogens, so about 57% are positive even
  at zero promotion hits, and about half of them are genotoxic (Ames-positive)
  initiators the promotion panel is not built to detect. Restricting to
  non-genotoxic (Ames-negative) carcinogens, the index still does not separate
  overall (AUC 0.51). The label rate by promotion-hit count is noisy and non-monotone
  on these small bands (about 52% at zero hits, 35% at one, 42% at two, 74% at three,
  67% at four or more), so any upward trend among the most active chemicals is
  suggestive at best. The lesson is anchor validity, not size.
- **The panel misses chronic inflammation (ext 3).** Mapped onto the ten key
  characteristics of carcinogens, the 12-assay panel reads receptor effects,
  oxidative stress, altered proliferation and (through the held-out DNA-damage
  assays) genotoxicity, so it covers 4 of 10 characteristics and has no readout for
  the other 6, including chronic inflammation, which is central to promotion.
  Measures of the same key characteristic do not agree more than measures of
  different ones (phi 0.17 within versus 0.20 across), so the textbook ontology
  does not predict co-activation either.
- **Connectivity nominates plausible promolytics (ext 5).** Querying L1000
  connectivity to reverse a promotion / inflammaging signature returns a shortlist
  led by an HSP90 inhibitor (geldanamycin), a Chk1 inhibitor (AZD-7762) and several
  EGFR inhibitors (canertinib, afatinib), with an HDAC inhibitor (vorinostat) and
  statins also recurring, all classes with chemoprevention literature. The mirror
  query for mimickers returns phorbol esters (PMA, ingenol dibenzoate) among the top
  hits, behind an unidentified Broad compound and narciclasine. The phorbol esters
  are the textbook tumour promoters, so their appearance is the positive control
  that the signature captures real promotion biology.

## Limits

- Non-genotoxic carcinogen is a proxy for tumour promoter, not an identity.
- The binary index measures breadth of pathway activation, not strength. Ext 1
  addresses this with qHTS potencies for the same pathways. ToxCast/invitroDB
  AC50 values would extend it further.
- The carcinogenicity anchor is either small (Lagunin) or large but selection-
  biased, with about half being genotoxic initiation the panel does not read
  (IARC/NTP). The closest-to-valid anchor is the non-genotoxic subset, and even
  there the index only separates outcomes at the top, not overall.
- The LINCS screen is hypothesis generation, not a claim about any single compound.

## Not yet built

**Toxicogenomic measures from Open TG-GATEs / DrugMatrix (ext 4).** This would add
expression-based promotion signatures as a second, orthogonal family of competing
measures, letting the audit test cross-platform agreement. It is deferred because
laptop-scale use of Open TG-GATEs means downloading and normalizing a large
processed expression archive (from the Life Science Database Archive at
dbarchive.biosciencedbc.jp), which is a project in itself. The LINCS connectivity
screen (ext 5) already brings in an orthogonal transcriptomic platform, so the
cross-platform check has a first foothold without it.
