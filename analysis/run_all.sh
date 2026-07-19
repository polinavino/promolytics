#!/usr/bin/env bash
# Run the whole promolytics pipeline in order. Every script re-teed its console
# output into outputs/<script>.txt, so after this finishes the tracked outputs/
# folder is the current result of the current code on the current data.
#
# Interpreter: the scripts need numpy, pandas, scipy, scikit-learn and rdkit.
# Set PYTHON to point at such an interpreter, or let this script find the local
# miniforge molml env that already has them.
set -euo pipefail
cd "$(dirname "$0")"

if [ -z "${PYTHON:-}" ]; then
  for cand in "$HOME/miniforge3/envs/molml/bin/python" \
              "$HOME/miniforge3/envs/neoantigen/bin/python" python3 python; do
    if command -v "$cand" >/dev/null 2>&1 && \
       "$cand" -c "import numpy,pandas,scipy,sklearn,rdkit" >/dev/null 2>&1; then
      PYTHON="$cand"; break
    fi
  done
fi
if [ -z "${PYTHON:-}" ]; then
  echo "No interpreter with the required packages found. Set PYTHON=..." >&2
  exit 1
fi
echo "using interpreter: $PYTHON"

echo "-- core audit --"
"$PYTHON" download_data.py
"$PYTHON" build_matrix.py
"$PYTHON" families.py
"$PYTHON" avg_rank.py
"$PYTHON" stability.py
"$PYTHON" anchor.py
"$PYTHON" honest_eval.py

# Extensions. The potency and anchor scripts make many PubChem calls and cache to
# data/raw, so they are slow on a cold cache and near-instant once cached. Run
# ext_potency before ext_anchor: they share the CID->InChIKey cache and must not
# run concurrently.
echo "-- extensions --"
"$PYTHON" ext_potency.py
"$PYTHON" ext_anchor.py
"$PYTHON" ext_kcc.py
"$PYTHON" ext_lincs.py
echo "done. see outputs/"
