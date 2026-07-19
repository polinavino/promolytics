"""Shared PubChem access helpers used by the extension scripts.

PubChem's public REST endpoints return 503 (ServerBusy) and 504 (Timeout)
freely under load, so every call here retries with exponential backoff, and the
CID to InChIKey mapping caches to disk incrementally so a re-run only fetches the
ids still missing. Both extension scripts share the same cache file.
"""

import io
import os
import time

import pandas as pd
import requests

import _common as C

IKEY_CACHE = os.path.join(C.RAW, "potency", "cid_inchikey.csv")


def get_with_retry(url, timeout=90, tries=5):
    """GET a URL, retrying transient PubChem errors. Returns a Response or None."""
    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code in (503, 504):
                raise requests.exceptions.HTTPError("transient %d" % r.status_code)
            r.raise_for_status()
            return r
        except (requests.exceptions.RequestException, requests.exceptions.HTTPError):
            time.sleep(2 ** (attempt + 1))
    return None


def _load_cache():
    if os.path.exists(IKEY_CACHE):
        c = pd.read_csv(IKEY_CACHE)
        return dict(zip(c["CID"].astype(int), c["ikey"].astype(str)))
    return {}


def _save_cache(cache):
    os.makedirs(os.path.dirname(IKEY_CACHE), exist_ok=True)
    pd.DataFrame({"CID": list(cache), "ikey": list(cache.values())}).to_csv(IKEY_CACHE, index=False)


def map_cid_to_ikey(cids, batch=100):
    """Map CIDs to InChIKey connectivity blocks, caching incrementally to disk."""
    cache = _load_cache()
    missing = [int(c) for c in cids if int(c) not in cache]
    failed = 0
    for i in range(0, len(missing), batch):
        chunk = missing[i:i + batch]
        ids = ",".join(str(c) for c in chunk)
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/InChIKey/CSV"
               % ids)
        r = get_with_retry(url)
        if r is None:
            failed += 1
            continue
        d = pd.read_csv(io.StringIO(r.text))
        for _, row in d.iterrows():
            cache[int(row["CID"])] = str(row["InChIKey"]).split("-")[0]
        if i % 1000 == 0:
            _save_cache(cache)
        time.sleep(0.3)
    _save_cache(cache)
    if failed:
        print("  %d CID batch(es) still unresolved, re-run to fill the gaps" % failed)
    return cache
