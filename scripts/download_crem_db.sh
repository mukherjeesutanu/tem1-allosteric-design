#!/usr/bin/env bash
# CReM fragment-replacement database (ChEMBL22, SA<=2.5, HAC<=12).
# ~846 MB compressed, 5.8 GB unpacked.
set -euo pipefail
DEST="$(dirname "$0")/../data/raw"
URL="https://zenodo.org/records/16909329/files/chembl22_sa25_hac12.db.gz?download=1"
mkdir -p "$DEST"
echo "downloading CReM database (~846 MB)…"
curl -L -C - -o "$DEST/chembl22_sa25_hac12.db.gz" "$URL"
gunzip -f "$DEST/chembl22_sa25_hac12.db.gz"
ls -lh "$DEST/chembl22_sa25_hac12.db"
