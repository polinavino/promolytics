#!/usr/bin/env python3
"""Extension 3: the panel seen through the key characteristics of carcinogens.

The predictive-tox and IARC community organises carcinogen mechanisms by the ten
key characteristics of carcinogens (KCC, Smith et al. 2016). Mapping the panel
onto that ontology does two things. It shows which key characteristics the Tox21
promotion panel actually covers, and it tests whether measures that share a key
characteristic agree more with each other than with measures of a different one,
which is what a well-formed mechanistic grouping should show.

The headline is a coverage gap: the panel reads receptor effects, oxidative
stress and altered proliferation, but has no readout for chronic inflammation or
immunosuppression, both central to promotion. That gap is why the connectivity
screen in ext_lincs (which reads the inflammatory program directly) is a genuine
complement, not a duplicate.
"""

import numpy as np

import _common as C
from families import phi_matrix

# Each assay mapped to its primary key characteristic of carcinogens.
ASSAY_KCC = {
    "NR-AhR": "KC8 receptor-mediated effects",
    "NR-PPAR-gamma": "KC8 receptor-mediated effects",
    "NR-ER": "KC8 receptor-mediated effects",
    "NR-ER-LBD": "KC8 receptor-mediated effects",
    "NR-AR": "KC8 receptor-mediated effects",
    "NR-AR-LBD": "KC8 receptor-mediated effects",
    "NR-Aromatase": "KC8 receptor-mediated effects",
    "SR-ARE": "KC5 oxidative stress",
    "SR-HSE": "KC5 oxidative stress",
    "SR-MMP": "KC10 altered cell proliferation/death",
    "SR-p53": "KC2 genotoxic",
    "SR-ATAD5": "KC2 genotoxic",
}

# The ten key characteristics, so we can name the ones with no panel readout.
ALL_KCC = [
    "KC1 electrophilic / metabolic activation",
    "KC2 genotoxic",
    "KC3 alters DNA repair / genomic instability",
    "KC4 epigenetic alterations",
    "KC5 oxidative stress",
    "KC6 chronic inflammation",
    "KC7 immunosuppression",
    "KC8 receptor-mediated effects",
    "KC9 immortalisation",
    "KC10 altered cell proliferation/death",
]


def main():
    C.rule("EXT 3: panel coverage of the key characteristics of carcinogens")
    m = C.load_matrix()

    # --- coverage ---
    covered = {}
    for a in C.ALL_ASSAYS:
        covered.setdefault(ASSAY_KCC[a], []).append(a)
    print("\nkey characteristics covered by the 12-assay panel:")
    for kcc in ALL_KCC:
        assays = covered.get(kcc, [])
        mark = "  " if assays else "**"
        print("  %s %-42s %s" % (mark, kcc, ", ".join(assays) if assays else "(no readout)"))
    gaps = [k for k in ALL_KCC if k not in covered]
    print("\ncoverage gaps (no assay reads these): %d of %d" % (len(gaps), len(ALL_KCC)))
    for g in gaps:
        print("  ** %s" % g)
    print("KC6 (chronic inflammation) is the most promotion-relevant gap.")

    # --- within-KCC vs cross-KCC agreement, on the promotion panel ---
    C.rule("Do measures of the same key characteristic agree more?")
    panel = C.PROMOTION_PANEL
    phi = phi_matrix(m, panel)
    within, cross = [], []
    for i in range(len(panel)):
        for j in range(i + 1, len(panel)):
            v = phi.iloc[i, j]
            if np.isnan(v):
                continue
            if ASSAY_KCC[panel[i]] == ASSAY_KCC[panel[j]]:
                within.append(v)
            else:
                cross.append(v)
    print("mean phi within the same key characteristic: %.3f (%d pairs)"
          % (np.mean(within), len(within)))
    print("mean phi across key characteristics:         %.3f (%d pairs)"
          % (np.mean(cross), len(cross)))
    print("=> measures of the same key characteristic agree %s than across"
          % ("more" if np.mean(within) > np.mean(cross) else "no more"))

    # --- per-chemical key-characteristic breadth ---
    C.rule("How many key characteristics does each chemical activate?")
    # a chemical activates a KCC if it hits any panel assay mapped to it
    kcc_of_panel = {}
    for a in panel:
        kcc_of_panel.setdefault(ASSAY_KCC[a], []).append(a)
    tested_all = m[panel].notna().all(axis=1)
    sub = m[tested_all]
    counts = np.zeros(len(sub), dtype=int)
    for kcc, assays in kcc_of_panel.items():
        hit_any = (sub[assays] == 1).any(axis=1).values
        counts += hit_any.astype(int)
    n_kcc = len(kcc_of_panel)
    print("complete-case chemicals: %d, spanning %d promotion key characteristics"
          % (len(sub), n_kcc))
    for k in range(n_kcc + 1):
        c = int((counts == k).sum())
        print("  activates %d of %d key characteristics: %5d chemicals (%.1f%%)"
              % (k, n_kcc, c, 100 * c / len(sub)))
    print("(a chemical hitting several key characteristics is a broader promoter")
    print(" than one hitting many assays inside a single characteristic)")


if __name__ == "__main__":
    with C.tee("ext_kcc"):
        main()
