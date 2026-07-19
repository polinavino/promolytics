#!/usr/bin/env python3
"""Extension 1: a potency-weighted promotability index from Tox21 qHTS AC50.

The binary index measures breadth of pathway activation. It cannot tell a weak
hit from a strong one, which is why so many chemicals end up as near-ties. This
extension replaces the binary hit-calls with AC50 potencies from the Tox21
quantitative high-throughput screen (qHTS), pulled from PubChem BioAssay, and
asks whether grading the measures by potency changes the ranking and breaks the
ties.

For each panel pathway we use the PubChem qHTS activation assay. Potency is read
as AC50 in micromolar for active substances and converted to pAC50 = 6 -
log10(AC50_uM), so a larger pAC50 is a more potent activator. Inactive-but-tested
substances get the bottom of that assay's scale, and untested substances stay
missing. The potency index is then the same average-percentile-rank aggregate as
the binary one, so the two are directly comparable.

Each column is labelled with the assay title PubChem actually returned, so the
names in the output are the names of the data, not an assumed mapping.
"""

import os
import time

import numpy as np
import pandas as pd

import _common as C
import pubchem_util as P

# PubChem qHTS activation assays for the promotion panel. Direction is activation
# of the pathway (agonist), except aromatase which PubChem screens as inhibition.
# The returned assay title is printed for verification so a wrong id is visible.
POTENCY_AIDS = {
    "NR-AhR": 743122,
    "NR-PPAR-gamma": 743140,
    "NR-ER": 743079,
    "NR-ER-LBD": 743077,
    "NR-AR": 743053,
    "NR-AR-LBD": 743040,
    "NR-Aromatase": 743139,   # aromatase inhibition (noted in output)
    "SR-ARE": 743219,
    "SR-HSE": 743228,
    "SR-MMP": 720637,
}

POT_RAW = os.path.join(C.RAW, "potency")
POT_MATRIX = os.path.join(C.DATA, "potency_matrix.csv")


def download_concise(aid):
    os.makedirs(POT_RAW, exist_ok=True)
    path = os.path.join(POT_RAW, "aid_%d.csv" % aid)
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%d/concise/CSV" % aid
    r = P.get_with_retry(url, timeout=120)
    if r is None or len(r.content) < 1000:
        raise RuntimeError("could not download AID %d after retries" % aid)
    with open(path, "wb") as fh:
        fh.write(r.content)
    time.sleep(0.5)
    return path


def parse_assay(aid):
    """Return (title, DataFrame[cid, pac50]) of active substances for one assay."""
    path = download_concise(aid)
    df = pd.read_csv(path, low_memory=False)
    title = df["Assay Name"].dropna().iloc[0] if "Assay Name" in df else ""
    active = df[df["Activity Outcome"] == "Active"].copy()
    active = active.dropna(subset=["CID", "Activity Value [uM]"])
    active["CID"] = active["CID"].astype(int)
    ac50 = pd.to_numeric(active["Activity Value [uM]"], errors="coerce")
    active = active.assign(ac50=ac50).dropna(subset=["ac50"])
    active = active[active["ac50"] > 0]
    active["pac50"] = 6.0 - np.log10(active["ac50"])
    # one value per compound (median across replicates)
    pac = active.groupby("CID")["pac50"].median()
    tested = set(df.dropna(subset=["CID"])["CID"].astype(int))
    return title, pac, tested


def percentile_index(mat, assays, min_tested=5):
    ranks = pd.DataFrame(index=mat.index)
    for a in assays:
        ranks[a] = mat[a].rank(method="average", pct=True)
    tested = mat[assays].notna().sum(axis=1)
    idx = ranks.mean(axis=1, skipna=True)
    idx[tested < min_tested] = np.nan
    return idx


def main():
    C.rule("EXT 1: potency-weighted promotability index (Tox21 qHTS AC50)")

    print("\nPulling qHTS AC50 data from PubChem (cached in data/raw/potency):")
    per_assay = {}
    tested_sets = {}
    titles = {}
    all_cids = set()
    for name, aid in POTENCY_AIDS.items():
        try:
            title, pac, tested = parse_assay(aid)
        except RuntimeError as e:
            print("  %-16s AID %-7d SKIPPED (%s)" % (name, aid, str(e)[:45]))
            continue
        per_assay[name] = pac
        tested_sets[name] = tested
        titles[name] = title
        all_cids |= set(pac.index) | tested
        print("  %-16s AID %-7d actives %5d  tested %5d" % (name, aid, len(pac), len(tested)))

    # a measure with no active substances (e.g. the aromatase inhibition summary,
    # which does not expose activation potency) carries no ranking signal, so drop it
    panel = [n for n in per_assay if len(per_assay[n]) > 0]
    dropped = [n for n in per_assay if len(per_assay[n]) == 0]
    if dropped:
        print("\ndropping %d measure(s) with no potency values: %s"
              % (len(dropped), ", ".join(dropped)))
    unavailable = [n for n in POTENCY_AIDS if n not in per_assay]
    if unavailable:
        print("note: %d measure(s) unavailable from PubChem this run: %s"
              % (len(unavailable), ", ".join(unavailable)))
    print("potency panel: %d measures" % len(panel))

    print("\nAssay titles returned by PubChem (verify against the column names):")
    for name in panel:
        print("  %-16s -> %s" % (name, titles[name][:66]))

    # map every CID to an InChIKey block
    print("\nmapping %d CIDs to InChIKey..." % len(all_cids))
    cid2ikey = P.map_cid_to_ikey(sorted(all_cids))

    # build wide potency matrix keyed on ikey: pAC50 for actives, 0.0 for
    # tested-inactive, NaN for untested
    ikeys = sorted({cid2ikey[c] for c in all_cids if c in cid2ikey})
    mat = pd.DataFrame(index=pd.Index(ikeys, name="ikey"), columns=panel, dtype=float)
    for name in panel:
        pac = per_assay[name]
        tested = tested_sets[name]
        # tested-inactive -> 0.0
        for c in tested:
            ik = cid2ikey.get(c)
            if ik is not None:
                mat.at[ik, name] = 0.0
        # actives -> pAC50 (overwrites the 0.0), median across CIDs sharing an ikey
        tmp = {}
        for cid, val in pac.items():
            ik = cid2ikey.get(cid)
            if ik is None:
                continue
            tmp.setdefault(ik, []).append(val)
        for ik, vals in tmp.items():
            mat.at[ik, name] = float(np.median(vals))

    mat = mat.reset_index()
    mat.to_csv(POT_MATRIX, index=False)
    print("wrote data/potency_matrix.csv (%d chemicals x %d potency measures)"
          % (len(mat), len(panel)))

    # potency index
    pot_idx = percentile_index(mat.set_index("ikey"), panel)
    pot_idx = pot_idx.dropna()
    print("\npotency index: %d chemicals scored (>= 5 measures tested)" % len(pot_idx))

    # compare to the binary index
    C.rule("Potency index vs binary index")
    binm = C.load_matrix()
    bin_idx = pd.read_csv(f"{C.DATA}/promotability_index.csv").set_index("ikey")["promotability"]
    common = pot_idx.index.intersection(bin_idx.index)
    from scipy.stats import spearmanr
    rho = spearmanr(pot_idx.loc[common], bin_idx.loc[common]).correlation
    print("chemicals in both indices: %d" % len(common))
    print("Spearman(potency, binary): %.3f" % rho)

    # how many binary ties does potency resolve?
    b = bin_idx.loc[common].round(4)
    tie_groups = b.groupby(b).groups
    tied_chems = sum(len(g) for g in tie_groups.values() if len(g) > 1)
    resolved = 0
    for val, grp in tie_groups.items():
        if len(grp) < 2:
            continue
        pv = pot_idx.loc[list(grp)]
        if pv.nunique() > 1:
            resolved += len(grp)
    print("chemicals sharing a binary-index value with another: %d" % tied_chems)
    print("of those, given a distinct value by potency:         %d (%.0f%%)"
          % (resolved, 100 * resolved / max(tied_chems, 1)))

    # top by potency
    C.rule("Top promotability by potency index")
    top = pot_idx.sort_values(ascending=False).head(10)
    matx = mat.set_index("ikey")
    for ik, v in top.items():
        active_paths = [a for a in panel if matx.at[ik, a] > 0]
        print("  %.3f  %s  active in %d pathways" % (v, ik, len(active_paths)))

    out = pot_idx.rename("potency_index").reset_index()
    out.to_csv(f"{C.DATA}/potency_index.csv", index=False)
    print("\nwrote data/potency_index.csv (%d chemicals)" % len(out))


if __name__ == "__main__":
    with C.tee("ext_potency"):
        main()
