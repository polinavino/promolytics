#!/usr/bin/env python3
"""Fetch the raw public datasets into data/raw/.

All three come from the Therapeutics Data Commons mirror on the Harvard
Dataverse, which is openly downloadable and needs no account. File ids are the
stable dataverse datafile ids used by PyTDC. Files already present are left
alone unless --force is passed.

  tox21                12 mechanistic hit-call assays + SMILES (competing measures)
  ames                 mutagenicity / genotoxicity label + SMILES (initiation proxy)
  carcinogens_lagunin  carcinogenicity label + SMILES (outcome anchor)
"""

import os
import sys

import requests

import _common as C

DATAVERSE = "https://dataverse.harvard.edu/api/access/datafile/"
FILES = {
    "tox21": 4259612,
    "ames": 4259564,
    "carcinogens_lagunin": 4259570,
}


def fetch(name, file_id, force=False):
    os.makedirs(C.RAW, exist_ok=True)
    dest = os.path.join(C.RAW, name + ".tab")
    if os.path.exists(dest) and not force:
        print("  %-22s present (%d bytes), skipping" % (name, os.path.getsize(dest)))
        return
    r = requests.get(DATAVERSE + str(file_id), timeout=180)
    r.raise_for_status()
    with open(dest, "wb") as fh:
        fh.write(r.content)
    print("  %-22s downloaded (%d bytes)" % (name, len(r.content)))


def main():
    force = "--force" in sys.argv
    C.rule("DOWNLOAD DATA: TDC datasets from Harvard Dataverse")
    for name, fid in FILES.items():
        fetch(name, fid, force=force)
    print("\nraw data in %s" % C.RAW.replace(C.ROOT + "/", ""))


if __name__ == "__main__":
    with C.tee("download_data"):
        main()
