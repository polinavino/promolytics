"""Shared helpers for the promolytics measurement-audit pipeline.

Two jobs:

1. Self-teeing output. Every analysis script wraps its body in `with tee(name):`
   so that everything printed to the console is also written, verbatim, to
   `outputs/<name>.txt`. That file is overwritten on every run, so the tracked
   copy in git is always the output of the current code on the current data.
   Outputs are the source of truth: keep them deterministic (no timestamps, no
   unseeded randomness) so a git diff of outputs/ shows real result changes only.

2. Paths and a loader for the harmonized matrix that build_matrix.py produces,
   so downstream scripts share one definition of the data.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "raw")
OUTPUTS = os.path.join(ROOT, "outputs")

# The harmonized chemical x assay matrix plus anchor labels, written by
# build_matrix.py and read by every downstream script.
MATRIX_CSV = os.path.join(DATA, "promolytics_matrix.csv")

# Fixed seed for every stochastic step in the pipeline. Deterministic outputs
# are the whole point of tracking outputs/ as source of truth.
SEED = 20260719

# The 12 Tox21 assay endpoints, grouped into mechanistic families. Each assay is
# a competing measure of promotional potential. The families are the pathway
# clusters the framework expects these measures to fall into.
#
# Selection rule (stated up front, see README): the promotability panel is the
# set of Tox21 endpoints that read promotion-relevant key characteristics of
# carcinogens: receptor-mediated signalling, oxidative/proteotoxic stress, and
# mitochondrial/proliferative stress. The two DNA-damage-response readouts
# (SR-p53, SR-ATAD5) are genotoxicity-leaning (initiation, not promotion) and are
# held out of the panel by default, but kept in the matrix so the audit can show
# they behave differently.
ASSAY_FAMILIES = {
    "nuclear_receptor": [
        "NR-AhR",         # aryl hydrocarbon receptor (xenobiotic, dioxin-like promotion)
        "NR-PPAR-gamma",  # peroxisome proliferator (classic non-genotoxic promotion)
        "NR-ER", "NR-ER-LBD",     # estrogen receptor (hormonal promotion)
        "NR-AR", "NR-AR-LBD",     # androgen receptor
        "NR-Aromatase",           # estrogen synthesis
    ],
    "oxidative_proteo_stress": [
        "SR-ARE",   # antioxidant response element (oxidative stress)
        "SR-HSE",   # heat-shock / proteotoxic stress
        "SR-MMP",   # mitochondrial membrane potential (proliferative/metabolic stress)
    ],
    "dna_damage": [   # genotoxicity-leaning, excluded from the promotability panel
        "SR-p53",
        "SR-ATAD5",
    ],
}

# The measures that make up the promotability panel (everything except the
# genotoxicity-leaning DNA-damage readouts).
PROMOTION_PANEL = ASSAY_FAMILIES["nuclear_receptor"] + ASSAY_FAMILIES["oxidative_proteo_stress"]
ALL_ASSAYS = PROMOTION_PANEL + ASSAY_FAMILIES["dna_damage"]

# Map each panel assay back to its family label.
ASSAY_TO_FAMILY = {}
for _fam, _assays in ASSAY_FAMILIES.items():
    for _a in _assays:
        ASSAY_TO_FAMILY[_a] = _fam


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


class tee:
    """Context manager: duplicate stdout to outputs/<name>.txt (overwritten)."""

    def __init__(self, name):
        os.makedirs(OUTPUTS, exist_ok=True)
        self.path = os.path.join(OUTPUTS, name + ".txt")
        self._fh = None
        self._old = None

    def __enter__(self):
        self._fh = open(self.path, "w")
        self._old = sys.stdout
        sys.stdout = _Tee(self._old, self._fh)
        return self

    def __exit__(self, *exc):
        sys.stdout = self._old
        self._fh.close()
        # Report where the output landed on the real console only.
        print("[written to %s]" % os.path.relpath(self.path, ROOT))
        return False


def rule(title=""):
    """A section divider, kept simple so outputs stay diff-friendly."""
    if title:
        print("\n" + "=" * 70)
        print(title)
        print("=" * 70)
    else:
        print("-" * 70)


def load_matrix():
    """Load the harmonized matrix, or fail with a clear message if not built yet."""
    import pandas as pd

    if not os.path.exists(MATRIX_CSV):
        raise SystemExit(
            "Harmonized matrix not found at %s.\nRun `python analysis/build_matrix.py` first."
            % MATRIX_CSV
        )
    return pd.read_csv(MATRIX_CSV)
