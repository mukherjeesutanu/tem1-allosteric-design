#!/usr/bin/env bash
# TEM-1 beta-lactamase structures used by the project.
#   1PZO  allosteric complex with CBT   (Horn & Shoichet, JMB 2004)
#   1PZP  allosteric complex with FTA
#   1BTL  apo TEM-1 (pocket-collapse control)
set -euo pipefail
DEST="$(dirname "$0")/../data/raw"
mkdir -p "$DEST"
for pdb in 1PZO 1PZP 1BTL; do
  curl -sL -o "$DEST/$pdb.pdb" "https://files.rcsb.org/download/$pdb.pdb"
  echo "$pdb: $(grep -c '^ATOM' "$DEST/$pdb.pdb") protein atoms"
done
