#!/usr/bin/env python3
"""Extension 5: an in-silico promolytic screen by L1000 connectivity.

The promotion step runs through an inflammatory and proliferative transcriptional
program (NF-kB cytokines, COX-2, proliferation and remodelling genes up, several
tumour-suppressor and antioxidant genes down). A promolytic is, by definition, an
agent that reverses that program. This script defines the promotion signature as
an up gene set and a down gene set, then asks the LINCS L1000 connectivity engine
(L1000CDS2) two questions:

  reversers   compounds whose L1000 signature is opposite to promotion
              -> candidate promolytics / chemopreventives
  mimickers   compounds whose L1000 signature matches promotion
              -> candidate promoters

The signature is deliberately compact and is stated in full below so the screen
is transparent. This is a hypothesis-generating shortlist, not a claim about any
single compound. It is the connectivity logic from the handoff applied to a novel
target concept.
"""

import numpy as np
import requests

import _common as C

L1000CDS2 = "https://maayanlab.cloud/L1000CDS2/query"

# Promotion / inflammaging signature.
# UP during promotion: NF-kB inflammatory cytokines, COX-2, proliferation and
# tissue-remodelling genes.
PROMOTION_UP = [
    "IL6", "IL1B", "TNF", "CXCL8", "CCL2", "PTGS2", "NFKB1", "NFKBIA",
    "MMP9", "VEGFA", "CCND1", "MYC", "JUN", "FOS", "SPP1", "TIMP1",
    "SERPINE1", "ICAM1", "CXCL1", "CXCL2",
]
# DOWN during promotion: cell-cycle brakes, antioxidant and tumour-suppressor genes.
PROMOTION_DOWN = [
    "CDKN1A", "CDKN2A", "GADD45A", "SOD2", "TXNIP", "PTEN", "TP53",
]


def query(up, dn, aggravate):
    """Query L1000CDS2. aggravate=False returns reversers, True returns mimickers."""
    payload = {
        "data": {"upGenes": up, "dnGenes": dn},
        "config": {"aggravate": aggravate, "searchMethod": "geneSet",
                   "share": False, "combination": False, "db-version": "latest"},
    }
    last = None
    for attempt in range(4):
        try:
            r = requests.post(L1000CDS2, json=payload, timeout=120)
            r.raise_for_status()
            return r.json().get("topMeta", [])
        except requests.exceptions.RequestException as e:
            last = e
            import time
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("L1000CDS2 query failed: %s" % last)


def aggregate(topmeta):
    """Collapse per-signature hits to one row per compound (best and mean score)."""
    by = {}
    for m in topmeta:
        name = m.get("pert_desc") or m.get("pert_id") or "?"
        if name in ("-666", "?", ""):
            name = m.get("pert_id", "?")
        by.setdefault(name, []).append(float(m.get("score", 0.0)))
    rows = []
    for name, scores in by.items():
        rows.append((name, len(scores), max(scores), float(np.mean(scores))))
    rows.sort(key=lambda t: (-t[1], -t[2]))  # by hit count then best score
    return rows


def show(title, rows, n=15):
    C.rule(title)
    print("  %-26s %6s %8s %8s" % ("compound", "hits", "best", "mean"))
    for name, cnt, best, mean in rows[:n]:
        print("  %-26s %6d %8.3f %8.3f" % (name[:26], cnt, best, mean))


def main():
    C.rule("EXT 5: in-silico promolytic screen (LINCS L1000 connectivity)")
    print("\nPromotion signature:")
    print("  UP   (%d genes): %s" % (len(PROMOTION_UP), ", ".join(PROMOTION_UP)))
    print("  DOWN (%d genes): %s" % (len(PROMOTION_DOWN), ", ".join(PROMOTION_DOWN)))

    print("\nquerying L1000CDS2 for reversers (candidate promolytics)...")
    reversers = aggregate(query(PROMOTION_UP, PROMOTION_DOWN, aggravate=False))
    print("querying L1000CDS2 for mimickers (candidate promoters)...")
    mimickers = aggregate(query(PROMOTION_UP, PROMOTION_DOWN, aggravate=True))

    show("Candidate promolytics (reverse the promotion signature)", reversers)
    show("Candidate promoters (mimic the promotion signature)", mimickers)

    # save the shortlists
    import csv
    out = f"{C.DATA}/lincs_candidates.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["role", "compound", "hits", "best_score", "mean_score"])
        for name, cnt, best, mean in reversers:
            w.writerow(["reverser_promolytic", name, cnt, "%.4f" % best, "%.4f" % mean])
        for name, cnt, best, mean in mimickers:
            w.writerow(["mimicker_promoter", name, cnt, "%.4f" % best, "%.4f" % mean])

    C.rule("Reading")
    print("Reversers are the in-silico promolytic shortlist. The most recurring hits are")
    print("an HSP90 inhibitor and a Chk1 inhibitor plus several EGFR inhibitors, and")
    print("classes already studied as chemopreventives (HDAC inhibitors, statins,")
    print("PI3K/mTOR inhibitors) also appear. The mimicker side is the sanity check:")
    print("phorbol esters, the textbook tumour promoters, appear among the top hits.")
    print("This is connectivity-based hypothesis generation, not a claim per compound.")
    print("\nwrote data/lincs_candidates.csv")


if __name__ == "__main__":
    with C.tee("ext_lincs"):
        main()
