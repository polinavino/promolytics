# Promolytics — handoff report, opportunity assessment, and starter plan

**Author-facing scoping document.** Written to decide whether to open a real project on computational
"promolytics" (early cancer prevention via the tumor-promotion step), and if so, where to start. It
assumes the reader knows the author's measurement-theory framework (competing measures of one latent
concept, consensus poset, average-rank canonical aggregate, desiderata) but not the promolytics field.
Style follows the author's rules (no semicolons, no em-dash pairs).

---

## 1. Verdict

**Promising, and an unusually good fit for the author's framework.** Start with a measurement and
prediction contribution on public proxy data, not with de-novo compound design (which is premature). The
single best hook is that the field is openly asking for a *measure* it does not yet have (a "promotion
Ames test"), which is the exact question the framework is built to answer.

Confidence: the biology and the funding are real and current. The computational side is nascent, which is
both the opportunity (low-hanging fruit) and the risk (definitions are unsettled, so part of the work is
defining the target).

---

## 2. What promolytics is (grounded)

Classical carcinogenesis has three stages: initiation (a mutation, genotoxic), promotion (expansion and
transformation of initiated cells, usually non-genotoxic, driven by inflammation, proliferation, and
signalling), and progression. "Promolytics" is a proposed class of compounds that block or reverse the
promotion step, so that initiated cells never become cancer. It is a cancer-interception idea.

The main current driver is the **PROMINENT** team, a Cancer Grand Challenges consortium (~$25M, UCSF,
Luke Gilbert and colleagues). Two facts matter for us:
- A stated goal is a screening test for the *promotional* ability of carcinogens, explicitly analogous to
  the Ames test for mutagenicity. Call it a "promotion assay" or "promotion Ames test." It does not exist
  yet.
- Their method so far is mostly wet-lab: genome-scale CRISPRi/a and Perturb-seq screens of tumour
  promoters and anti-promoting factors. New findings were posted October 2025.

Implication: the in-silico promolytics space is close to empty. Adjacent computational fields
(non-genotoxic carcinogen prediction, toxicogenomics, chemoprevention repurposing) are where the usable
data and methods already live.

Caveat to keep visible: "non-genotoxic carcinogen" is broader than "tumour promoter." Promoters are a
subset. Public data proxies promotion imperfectly, and any analysis must say so.

---

## 3. Why the framework fits unusually well

The framework's home question is: several formally distinct measures each claim to quantify one latent
concept, they disagree, and the field benchmarks them to pick a winner instead of mapping where the
concept decides. Promolytics gives this question in an unusually pure and *open* form.

- The latent concept is **promotional potential** (promotability) of a chemical or a perturbation.
- The candidate measures already exist as **mechanistic assays**: NF-kB and inflammatory readouts,
  nuclear-receptor activation (AhR, CAR, PXR, ER), oxidative-stress and proliferation assays in
  Tox21/ToxCast, plus various proposed non-genotoxic-carcinogen scores and toxicogenomic biomarker panels.
- These measures will disagree, will cluster into pathway families, and no single one is canonical. That
  is the consensus-poset and average-rank setting exactly.
- Better still, the field wants a **desiderata answer**: what should a well-formed promotion assay
  satisfy? That is the author's D-axiom / G-axiom method applied to a live, unmet need. Unlike selectivity
  or clocks (where measures pre-exist and the contribution is to reconcile them), here the measure itself
  is being defined, so there is room to shape the scaffold rather than only audit it.

This is arguably a stronger fit than any of the six domains already in the meta-paper, because the target
field is explicitly stuck on a measurement question.

---

## 4. Public data landscape (concrete)

No dedicated "promolytics dataset" exists. The usable substrate is adjacent public data.

- **Tox21** (public, NIH/EPA/FDA/NCATS). ~8,000 chemicals across a panel of nuclear-receptor and
  stress-response assays. Bulk download and via the Tox21 data browser. Many endpoints proxy
  promotion-relevant pathways.
- **ToxCast / invitroDB** (US EPA). ~4,700 chemicals across hundreds of in-vitro assay endpoints
  (inflammation, receptors, proliferation, cytotoxicity). Public, downloadable as a database release.
- **NICEATM / NIEHS carcinogenicity resources.** A public effort has imputed Tox21/ToxCast results across
  ~10,000 chemicals and ~2,000 assays mapped to the "key characteristics of carcinogens" (KCC). This is
  the closest thing to a promotion-mechanism feature matrix.
- **Carcinogenicity labels.** CPDB (Carcinogenic Potency Database), IARC monograph classifications, and
  curated genotoxic vs non-genotoxic splits. The genotoxic/non-genotoxic label is the practical proxy for
  initiator vs promoter.
- **Genotoxicity / Ames** data (public, e.g. via the Therapeutics Data Commons and Hansen benchmark) to
  *subtract* the initiation signal and isolate promotion-like non-genotoxic carcinogens.
- **Toxicogenomics.** TG-GATEs (Open TG-GATEs, public, rat and human hepatocyte expression under many
  compounds and doses) and DrugMatrix. These give expression signatures of promotion-relevant exposures.
- **LINCS L1000 / Connectivity Map** (public). Perturbation transcriptomics for thousands of compounds,
  usable to search for compounds that reverse an inflammation or promotion signature (an in-silico
  promolytic screen).
- **Therapeutics Data Commons (TDC)** packages several of the tox datasets in ML-ready form, which lowers
  the start-up cost.

Access reality: Tox21, ToxCast/invitroDB, TG-GATEs, LINCS, and TDC are openly downloadable and
laptop-scale to medium-scale. CPDB/IARC labels need light curation. Nothing here is gated.

---

## 5. Formal-analysis opportunities

Ordered from most framework-native to most standard.

1. **A measurement audit of promotability proxies (framework-native).** Take the mechanistic assays that
   each proxy promotion, treat them as competing measures of one concept, and run the protocol: families
   by induced order, consensus poset, average-rank canonical "promotability index," near-tie analysis, and
   an external anchor (non-genotoxic carcinogenicity outcome). Deliverable: a principled, weight-free
   promotability score plus a map of where the assays agree and where the concept is undecided. Framed as
   "toward a well-formed promotion assay."
2. **Desiderata for a promotion assay (conceptual).** State what a well-formed promotion measure should
   satisfy: a reliability/domain gate, mechanism/type declaration, monotonicity against an
   outcome anchor, ordinal stability under nuisance, cross-platform reproducibility, and intervention
   consistency. This is cheap, distinctive, and directly useful to the PROMINENT community.
3. **Non-genotoxic carcinogen prediction with honest evaluation (ADMET-sibling).** Predict non-genotoxic
   carcinogenicity from mechanistic assay fingerprints, and show the honest-split (scaffold or
   chemical-cluster holdout) performance drop versus random cross-validation. This reuses the author's
   evaluation-overstatement result in a new, high-stakes setting.
4. **In-silico promolytic screen via connectivity (LINCS).** Define a promotion or inflammaging
   transcriptomic signature, then rank compounds whose LINCS signature reverses it, to nominate candidate
   promolytics or chemopreventives. Standard connectivity logic, novel target concept.

---

## 6. Low-hanging fruit (prioritized, with feasibility)

1. **Promotability measurement audit (HIGH fit, MEDIUM effort).** Everything needed is public
   (Tox21/ToxCast assays + non-genotoxic carcinogen labels). Reuses the existing `avg_rank.py`,
   families, near-tie, and anchor machinery almost unchanged. This is the flagship first project and the
   one that plants a flag on the field's open measurement question.
2. **Desiderata note (HIGH fit, LOW effort).** A short, opinionated document proposing what a promotion
   assay must satisfy, with each axiom tied to a concrete failure mode seen in the audit. Pairs naturally
   with project 1 and is a fast way to be visible to PROMINENT.
3. **NGC honest-evaluation benchmark (MEDIUM fit, MEDIUM effort).** A clean, reusable benchmark showing
   that promotion-relevant carcinogen predictors are overstated under random CV. Useful and citable on its
   own.
4. **LINCS promolytic screen (MEDIUM fit, MEDIUM-HIGH effort).** Heavier data handling, but a memorable
   "candidate promolytics" shortlist is a good communication object.

Recommended sequence: 1 and 2 together as one repo and one short paper, then 3 as a companion, then 4 if
momentum is there.

---

## 7. Honest caveats and risks

- The field is embryonic and consortium-led. Definitions of a promotion measure are unsettled, so part of
  the contribution is defining the target. That is opportunity and exposure at once.
- Non-genotoxic carcinogen is a proxy for tumour promoter, not an identity. State this and, where
  possible, restrict to agents with promotion evidence rather than all non-genotoxic carcinogens.
- Outcome labels are noisy and sparse (animal bioassays, IARC categories). The anchor is crude, which the
  framework tolerates but the write-up must flag.
- A pure in-silico promolytic *drug* is out of near-term reach. The tractable contribution is measurement,
  prediction, and prioritization, not de-novo molecule design.
- Standing on PROMINENT's shoulders: engage the community rather than appear to scoop it. A measurement
  scaffold plus an honest benchmark is complementary to their wet-lab program, not competitive.

---

## 8. Suggested first project and repo scaffold

If pursued, mirror the meta-repo conventions (per-analysis scripts, self-teeing `outputs/`, stored numbers
as source of truth, a README documenting data provenance and the measure-selection rule).

```
promolytics/
  README.md                     # what/why, data provenance, selection rule, results, limits
  data/                         # Tox21/ToxCast assay matrix, NGC labels, curated promoter set
  analysis/
    build_matrix.py             # chemical x assay matrix + non-genotoxic-carcinogen labels
    families.py                 # competing promotability measures: correlation, families
    avg_rank.py                 # consensus poset + average-rank promotability index (reuse meta-repo)
    anchor.py                   # orientation + anchor (carcinogenicity outcome), extremes vs middle
    honest_eval.py              # NGC predictor: random CV vs scaffold/cluster holdout
    outputs/                    # tracked console outputs (source of truth)
  DESIDERATA.md                 # what a well-formed promotion assay should satisfy
  HANDOFF.md                    # this document
```

First concrete steps:
1. Acquire Tox21 + ToxCast/invitroDB assay matrices and a non-genotoxic vs genotoxic carcinogen label set
   (CPDB/IARC + an Ames source), harmonized on a chemical key (InChIKey or CASRN).
2. Select the promotion-relevant assays (inflammation, nuclear receptors, proliferation, oxidative stress)
   as the competing-measures panel, stating the selection rule up front.
3. Run the measurement audit (families, consensus, average-rank, near-tie, anchor) and draft DESIDERATA.md
   from what the audit reveals.

---

## 9. People and venues (for later, not now)

- PROMINENT / Cancer Grand Challenges (Luke Gilbert, UCSF and collaborators) for the biology and the
  promotion-assay goal.
- The predictive-tox and IARC key-characteristics-of-carcinogens community (NICEATM/NIEHS, US EPA CompTox)
  for the assay and carcinogenicity data and for an audience that already argues about mechanistic scores.
- Venues that fit a measurement or benchmark paper here: a computational-toxicology or cheminformatics
  journal for projects 1 and 3, and the author's usual measurement-theory framing for the desiderata.

---

## 10. Sources
- PROMINENT team, Cancer Grand Challenges: https://www.cancergrandchallenges.org/prominent
- Luke Gilbert lab, UCSF: https://cancer.ucsf.edu/people/gilbert.luke
- Non-genotoxic carcinogen consensus biomarkers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5259716/
- Toxicophore-informed ML on Tox21 for organ-specific carcinogenicity:
  https://www.sciencedirect.com/science/article/abs/pii/S0269749125018937
- NIEHS/NICEATM AI carcinogenicity models (Tox21/ToxCast imputation):
  https://ntp.niehs.nih.gov/iccvamreport/2023/technology/comp-tools-dev/22-niehs-carc-models
- Therapeutics Data Commons toxicity tasks: https://tdcommons.ai/single_pred_tasks/tox/
- Proactive prevention and cancer interception review:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12485381/
