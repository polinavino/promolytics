#!/usr/bin/env python3
"""Honest evaluation: random cross-validation overstates structure-based models.

A recurring result in the author's framework is that predictive-tox models look
better under random splits than under splits that hold out whole chemical
scaffolds. Random folds leak near-identical analogues across the train/test line,
so the reported number rewards memorising scaffolds rather than generalising to
new chemistry. This script reproduces that gap in the promotion setting.

For each target we predict the label from an ECFP4 structural fingerprint with a
random forest, and compare two evaluations:

  random    5-fold stratified cross-validation (the usual, optimistic number)
  scaffold  5-fold grouped by Bemis-Murcko scaffold so no scaffold is split

Targets:
  NR-AhR      a well-powered promotion-relevant mechanism (aryl hydrocarbon receptor)
  SR-ARE      a well-powered promotion-relevant mechanism (oxidative stress)
  carcinogen  the outcome anchor itself (small, from the Lagunin set)

The honest number is the scaffold one. The gap is the overstatement.
"""

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import roc_auc_score

import _common as C

RDLogger.DisableLog("rdApp.*")

N_BITS = 2048
RADIUS = 2
N_TREES = 200


def fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=np.int8)
    from rdkit.DataStructs import ConvertToNumpyArray
    ConvertToNumpyArray(fp, arr)
    return arr


def scaffold_of(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return None
    return scaf if scaf else "acyclic"


def cv_auc(X, y, splits):
    """Mean and sd of held-out AUC over the given (train, test) index splits."""
    aucs = []
    for tr, te in splits:
        if len(np.unique(y[te])) < 2 or len(np.unique(y[tr])) < 2:
            continue
        clf = RandomForestClassifier(
            n_estimators=N_TREES, n_jobs=-1, random_state=C.SEED, class_weight="balanced"
        )
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
    return float(np.mean(aucs)), float(np.std(aucs)), len(aucs)


def evaluate(name, smiles, labels):
    """Featurize, then compare random vs scaffold cross-validation."""
    df = pd.DataFrame({"smiles": smiles, "y": labels}).dropna()
    df = df[df["y"].isin([0, 1])]
    fps, scafs, ys = [], [], []
    for s, y in zip(df["smiles"], df["y"]):
        fp = fingerprint(s)
        if fp is None:
            continue
        fps.append(fp)
        scafs.append(scaffold_of(s) or "acyclic")
        ys.append(int(y))
    X = np.vstack(fps)
    y = np.array(ys)
    scafs = np.array(scafs)

    n_pos = int(y.sum())
    n_scaf = len(set(scafs))
    print("\n%s: %d compounds (%d positive, %.1f%%), %d distinct scaffolds"
          % (name, len(y), n_pos, 100 * n_pos / len(y), n_scaf))

    # random stratified 5-fold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=C.SEED)
    rand_splits = list(skf.split(X, y))
    r_mean, r_sd, r_k = cv_auc(X, y, rand_splits)

    # scaffold-grouped 5-fold (no scaffold spans the split)
    gkf = GroupKFold(n_splits=5)
    scaf_splits = list(gkf.split(X, y, groups=scafs))
    s_mean, s_sd, s_k = cv_auc(X, y, scaf_splits)

    print("  random    AUC %.3f +/- %.3f  (%d folds)" % (r_mean, r_sd, r_k))
    print("  scaffold  AUC %.3f +/- %.3f  (%d folds)" % (s_mean, s_sd, s_k))
    print("  overstatement (random - scaffold): %+.3f" % (r_mean - s_mean))
    return r_mean, s_mean


def main():
    C.rule("HONEST EVAL: random vs scaffold-held-out cross-validation")
    m = C.load_matrix()

    results = []
    for target in ("NR-AhR", "SR-ARE"):
        r, s = evaluate(target, m["smiles"], m[target])
        results.append((target, r, s))

    # carcinogen: use the full Lagunin set (more positives than the Tox21 overlap)
    carc = pd.read_csv(f"{C.RAW}/carcinogens_lagunin.tab", sep="\t")
    r, s = evaluate("carcinogen", carc["Drug"], carc["Y"])
    results.append(("carcinogen", r, s))

    C.rule("Summary")
    print("  %-12s %8s %8s %12s" % ("target", "random", "scaffold", "overstated"))
    for name, r, s in results:
        print("  %-12s %8.3f %8.3f %+12.3f" % (name, r, s, r - s))
    print("\nEvery target loses accuracy under the honest scaffold split. A promotion")
    print("model reported on random folds overstates how well it will call the")
    print("chemistry it has not seen. Report the scaffold number.")


if __name__ == "__main__":
    with C.tee("honest_eval"):
        main()
