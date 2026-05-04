#!/bin/bash
set -euo pipefail

BASE="$(cd "$(dirname "$0")" && pwd)"
DATASET="${1:-$BASE/tests/al_qha_el_litfix_new}"
OUT_SUMMARY="$DATASET/qha_elastic_summary.drop_v04.in"

if [ ! -f "$DATASET/qha_elastic_summary.in" ]; then
    echo "Missing summary file: $DATASET/qha_elastic_summary.in" >&2
    exit 1
fi

python3 - "$DATASET" "$OUT_SUMMARY" <<'PY'
import sys
from pathlib import Path

dataset = Path(sys.argv[1])
out_summary = Path(sys.argv[2])
source = dataset / "qha_elastic_summary.in"

energies = {}
for volume in ["v01", "v02", "v03", "v05", "v06", "v07", "v08", "v09"]:
    qe_out = dataset / volume / "qe.out"
    energy = None
    for line in qe_out.read_text().splitlines():
        if "!" in line and "total energy" in line:
            parts = line.split()
            if len(parts) >= 5:
                energy = parts[4]
    if energy is None:
        raise SystemExit(f"Missing SCF energy for {volume}")
    energies[volume] = energy

output_lines = []
for line in source.read_text().splitlines():
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        output_lines.append(line)
        continue

    parts = stripped.split()
    if len(parts) < 5:
        output_lines.append(line)
        continue

    vol, _, scf_template, elastic_dir, dos_path = parts[:5]
    volume = scf_template.split("/")[0]
    if volume == "v04":
        continue
    energy = energies.get(volume, parts[1])
    output_lines.append(f"{vol:<17} {energy:<23} {scf_template:<20} {elastic_dir:<22} {dos_path}")

out_summary.write_text("\n".join(output_lines) + "\n")
PY

cd "$DATASET"
"$BASE/build/qepp" qha_elastic -post "$(basename "$OUT_SUMMARY")" --tmin 0 --tmax 1000 --dt 100
