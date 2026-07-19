#!/usr/bin/env python3
"""Build the harmonized chemical x assay matrix plus outcome anchors.

Inputs (data/raw/, downloaded from the TDC Harvard Dataverse mirror):
  tox21.tab               12 mechanistic hit-call assays + SMILES  (competing measures)
  ames.tab                mutagenicity / genotoxicity label + SMILES (initiation proxy)
  carcinogens_lagunin.tab carcinogenicity label + SMILES            (outcome anchor)

Output (data/promolytics_matrix.csv), one row per chemical keyed on the InChIKey
connectivity block (first 14 chars, salt-stripped), with:
  - the 12 Tox21 assay hit-calls (0/1/NaN)
  - ames        genotoxicity label (0/1/NaN)
  - carcinogen  carcinogenicity label (0/1/NaN)   <- outcome anchor
  - ngc         non-genotoxic carcinogen = carcinogen & not mutagenic (0/1/NaN)

The InChIKey connectivity block ignores stereochemistry and protonation, which is
the right granularity for joining assay records to bioassay labels across sources.
"""

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

import _common as C

RDLogger.DisableLog("rdApp.*")  # silence rdkit parse warnings, we count failures ourselves


def largest_fragment(mol):
    """Return the largest connected fragment (drops salts/counterions)."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        return mol
    return max(frags, key=lambda m: m.GetNumAtoms())


def inchikey_block(smiles):
    """SMILES -> InChIKey connectivity block (14 chars), or None on failure."""
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = largest_fragment(mol)
        key = Chem.MolToInchiKey(mol)
    except Exception:
        return None
    if not key:
        return None
    return key.split("-")[0]


def keyed(df, smiles_col):
    """Add an 'ikey' connectivity-block column, report and drop parse failures."""
    df = df.copy()
    df["ikey"] = df[smiles_col].map(inchikey_block)
    n_fail = df["ikey"].isna().sum()
    print("  %-24s %5d rows, %4d unparsed/failed keys" % (smiles_col, len(df), n_fail))
    return df.dropna(subset=["ikey"])


def collapse(df, value_cols, how="max"):
    """Collapse duplicate InChIKeys to one row per chemical.

    For binary assay/label columns we take the max (a chemical counts as a hit /
    positive if any record for it is a hit / positive). This is the standard
    permissive collapse for Tox21-style replicate calls.
    """
    agg = {c: how for c in value_cols}
    return df.groupby("ikey", as_index=False).agg(agg)


def main():
    C.rule("BUILD MATRIX: harmonize Tox21 assays + genotoxicity + carcinogenicity")

    tox = pd.read_csv(f"{C.RAW}/tox21.tab", sep="\t")
    ames = pd.read_csv(f"{C.RAW}/ames.tab", sep="\t")
    carc = pd.read_csv(f"{C.RAW}/carcinogens_lagunin.tab", sep="\t")

    print("\nRaw sources:")
    print("  tox21        %5d rows x %d assays" % (len(tox), len(C.ALL_ASSAYS)))
    print("  ames         %5d rows (%d mutagenic / %d not)"
          % (len(ames), int((ames.Y == 1).sum()), int((ames.Y == 0).sum())))
    print("  carcinogens  %5d rows (%d carcinogen / %d not)"
          % (len(carc), int((carc.Y == 1).sum()), int((carc.Y == 0).sum())))

    print("\nAssigning InChIKey connectivity blocks:")
    tox = keyed(tox, "X")
    ames = keyed(ames, "Drug")
    carc = keyed(carc, "Drug")

    # Collapse each source to one row per chemical.
    tox_keep = C.ALL_ASSAYS
    tox_c = collapse(tox, tox_keep, how="max")
    # keep a representative SMILES per key for later scaffold splitting
    smiles_map = tox.drop_duplicates("ikey").set_index("ikey")["X"]
    tox_c["smiles"] = tox_c["ikey"].map(smiles_map)

    ames_c = collapse(ames.rename(columns={"Y": "ames"}), ["ames"], how="max")
    carc_c = collapse(carc.rename(columns={"Y": "carcinogen"}), ["carcinogen"], how="max")

    print("\nUnique chemicals after collapse:")
    print("  tox21 %d, ames %d, carcinogens %d"
          % (len(tox_c), len(ames_c), len(carc_c)))

    # Left-join labels onto the assay matrix. The assay panel is the substrate of
    # the audit, so the matrix is anchored on Tox21 chemicals.
    m = tox_c.merge(ames_c, on="ikey", how="left").merge(carc_c, on="ikey", how="left")

    # Non-genotoxic carcinogen: carcinogenic AND not mutagenic. Defined only where
    # both labels are known.
    both_known = m["carcinogen"].notna() & m["ames"].notna()
    ngc = np.where(both_known, ((m["carcinogen"] == 1) & (m["ames"] == 0)).astype(float), np.nan)
    m["ngc"] = ngc

    # Order columns.
    cols = ["ikey", "smiles"] + C.ALL_ASSAYS + ["ames", "carcinogen", "ngc"]
    m = m[cols]
    m.to_csv(C.MATRIX_CSV, index=False)

    C.rule("Harmonized matrix summary")
    print("chemicals (Tox21-anchored): %d" % len(m))
    print("with any assay call:        %d" % (m[C.ALL_ASSAYS].notna().any(axis=1).sum()))
    print("with ames label:            %d" % m["ames"].notna().sum())
    print("with carcinogen label:      %d" % m["carcinogen"].notna().sum())
    print("  of which carcinogenic:    %d" % int((m["carcinogen"] == 1).sum()))
    print("with ngc defined:           %d" % m["ngc"].notna().sum())
    print("  of which non-genotoxic carcinogens: %d" % int((m["ngc"] == 1).sum()))

    print("\nPromotability panel (competing measures), hit rate among tested:")
    for a in C.PROMOTION_PANEL:
        col = m[a]
        print("  %-16s tested %5d  hits %4d  (%.1f%%)  [%s]"
              % (a, col.notna().sum(), int((col == 1).sum()),
                 100 * (col == 1).sum() / max(col.notna().sum(), 1),
                 C.ASSAY_TO_FAMILY[a]))
    print("\nHeld out (genotoxicity-leaning, not in panel):")
    for a in C.ASSAY_FAMILIES["dna_damage"]:
        col = m[a]
        print("  %-16s tested %5d  hits %4d  (%.1f%%)"
              % (a, col.notna().sum(), int((col == 1).sum()),
                 100 * (col == 1).sum() / max(col.notna().sum(), 1)))

    print("\nwrote %s" % C.MATRIX_CSV.replace(C.ROOT + "/", ""))


if __name__ == "__main__":
    with C.tee("build_matrix"):
        main()
