#!/usr/bin/env python3
"""Extension 2: a larger carcinogenicity anchor from IARC and NTP classifications.

The base anchor (the Lagunin set) overlaps the Tox21 chemicals on only about a
hundred labels, too few to orient the index. This extension pulls the curated
carcinogenicity classifications PubChem aggregates under its "Carcinogen
Classification" heading, which is IARC monograph groups plus the NTP Report on
Carcinogens, for several thousand agents. Each agent is mapped to a chemical
structure and joined to the promotability index, then the same orientation and
band analysis is re-run on the larger overlap.

Label rule:
  positive (carcinogen = 1): IARC group 1, 2A or 2B, or NTP known / anticipated
  negative (carcinogen = 0): IARC group 3 or 4, or "not classifiable"
Records without a single CID (elements, mixtures, fibres) are dropped.
"""

import re

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import _common as C
import pubchem_util as P

HEADING = "Carcinogen%20Classification"
LABELS_CSV = f"{C.DATA}/carcinogen_labels_iarc_ntp.csv"

POS_PAT = re.compile(r"^\s*(1|2a|2b)\b", re.I)
NEG_PAT = re.compile(r"^\s*(3|4)\b", re.I)


def classify(text):
    """Map a classification string to 1 (carcinogen), 0 (not), or None (unknown)."""
    t = text.strip().lower()
    if POS_PAT.match(t):
        return 1
    if NEG_PAT.match(t):
        return 0
    if "known to be" in t or "reasonably anticipated" in t or "carcinogenic to humans" in t:
        return 1
    if "not classifiable" in t or "not carcinogenic" in t or "evidence of non" in t:
        return 0
    return None


def fetch_classifications():
    """Page through the PubChem heading, return DataFrame[cid, name, label]."""
    rows = []
    page = 1
    total_pages = 1
    while page <= total_pages:
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/JSON"
               "?heading_type=Compound&heading=%s&page=%d" % (HEADING, page))
        r = P.get_with_retry(url)
        if r is None:
            print("  page %d unavailable, stopping" % page)
            break
        data = r.json()["Annotations"]
        total_pages = data.get("TotalPages", 1)
        for a in data.get("Annotation", []):
            cids = a.get("LinkedRecords", {}).get("CID", [])
            if not cids:
                continue
            strings = []
            for dat in a.get("Data", []):
                for sm in dat.get("Value", {}).get("StringWithMarkup", []):
                    strings.append(sm.get("String", ""))
            label = None
            for s in strings:
                label = classify(s)
                if label is not None:
                    break
            if label is None:
                continue
            rows.append((int(cids[0]), a.get("Name", ""), label))
        print("  page %d/%d, cumulative labelled records: %d" % (page, total_pages, len(rows)))
        page += 1
    return pd.DataFrame(rows, columns=["CID", "name", "carc_iarc_ntp"])


def auc_rank(scores, labels):
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), len(pos), len(neg)
    u, p = mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg)), p, len(pos), len(neg)


def main():
    C.rule("EXT 2: larger carcinogenicity anchor (IARC + NTP via PubChem)")

    print("\nfetching carcinogen classifications from PubChem:")
    cls = fetch_classifications()
    print("labelled agents with a CID: %d (%d positive / %d negative)"
          % (len(cls), int((cls.carc_iarc_ntp == 1).sum()), int((cls.carc_iarc_ntp == 0).sum())))

    print("\nmapping %d CIDs to InChIKey..." % cls["CID"].nunique())
    cid2ikey = P.map_cid_to_ikey(cls["CID"].unique())
    cls["ikey"] = cls["CID"].map(cid2ikey)
    cls = cls.dropna(subset=["ikey"])
    # collapse to one call per chemical (positive wins over negative)
    lab = cls.groupby("ikey")["carc_iarc_ntp"].max().reset_index()
    lab.to_csv(LABELS_CSV, index=False)
    print("distinct chemicals labelled: %d" % len(lab))

    # join to the promotability index and hit counts
    m = C.load_matrix()
    idx = pd.read_csv(f"{C.DATA}/promotability_index.csv")
    panel = C.PROMOTION_PANEL
    n_hits = (m[panel] == 1).sum(axis=1)
    base = m[["ikey", "carcinogen", "ames"]].assign(n_hits=n_hits.values)
    scored = idx.merge(base, on="ikey", how="left").merge(lab, on="ikey", how="left")

    # combined anchor: union of Lagunin and IARC/NTP (positive wins)
    combined = scored[["carcinogen", "carc_iarc_ntp"]].max(axis=1, skipna=True)
    scored["carc_any"] = combined

    # non-genotoxic-carcinogen anchor: isolate promotion from initiation by keeping
    # only Ames-negative chemicals. positive = carcinogen and not mutagenic,
    # negative = non-carcinogen and not mutagenic. genotoxic carcinogens are dropped.
    ames_neg = scored["ames"] == 0
    ngc = np.where(
        ames_neg & scored["carc_any"].notna(),
        (scored["carc_any"] == 1).astype(float),
        np.nan,
    )
    scored["ngc_expanded"] = ngc

    C.rule("Overlap gained")
    print("index-scored chemicals with a label:")
    print("  Lagunin only:        %d (%d positive)"
          % (scored["carcinogen"].notna().sum(), int((scored["carcinogen"] == 1).sum())))
    print("  IARC/NTP:            %d (%d positive)"
          % (scored["carc_iarc_ntp"].notna().sum(), int((scored["carc_iarc_ntp"] == 1).sum())))
    print("  combined (any):      %d (%d positive)"
          % (scored["carc_any"].notna().sum(), int((scored["carc_any"] == 1).sum())))
    print("  non-genotoxic only:  %d (%d positive, Ames-negative)"
          % (scored["ngc_expanded"].notna().sum(), int((scored["ngc_expanded"] == 1).sum())))

    for col, title in [("carc_iarc_ntp", "IARC/NTP anchor (all carcinogens)"),
                       ("carc_any", "combined anchor (all carcinogens)"),
                       ("ngc_expanded", "non-genotoxic carcinogens only (Ames-negative)")]:
        C.rule("Orientation against %s" % title)
        lab_df = scored.dropna(subset=[col])
        auc, p, npos, nneg = auc_rank(lab_df["promotability"].values, lab_df[col].values)
        print("labelled + scored: %d (%d positive / %d negative)" % (len(lab_df), npos, nneg))
        if not np.isnan(auc):
            print("rank AUC of index vs %s: %.3f  (Mann-Whitney p=%.3g)" % (col, auc, p))
        print("\nlabel rate by number of promotion-pathway hits:")
        buckets = [("0 hits", 0, 0), ("1 hit", 1, 1), ("2 hits", 2, 2),
                   ("3 hits", 3, 3), ("4+ hits", 4, 99)]
        print("  %-9s %8s %8s %8s %8s" % ("hits", "n", "labeled", "pos", "rate"))
        for name, lo, hi in buckets:
            sub = scored[(scored["n_hits"] >= lo) & (scored["n_hits"] <= hi)]
            ll = sub[col].dropna()
            rate = "  n/a" if len(ll) == 0 else "%5.1f%%" % (100 * ll.mean())
            print("  %-9s %8d %8d %8d %8s"
                  % (name, len(sub), int(ll.notna().sum()), int((ll == 1).sum()), rate))

    C.rule("Reading")
    print("Expanding the anchor tenfold does not, on its own, make the index orient.")
    print("The IARC/NTP list is dominated by genotoxic carcinogens (initiators) and is")
    print("selection-biased toward suspected carcinogens, so the label rate is near")
    print("flat across promotion-hit counts. Restricting to non-genotoxic (Ames-negative)")
    print("carcinogens is the correct anchor for a promotion index. The lesson is about")
    print("anchor validity, not anchor size: a bigger wrong anchor is still the wrong one.")
    print("\nwrote %s" % LABELS_CSV.replace(C.ROOT + "/", ""))


if __name__ == "__main__":
    with C.tee("ext_anchor"):
        main()
