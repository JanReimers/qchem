#!/usr/bin/env python3
"""Extract the relativistic DHF data from DHF_GS_Energies.csv into JSON.

Keeps only Z, Conf, the relativistic total energy (TE_rel) and the
relativistic orbital eigenvalues (the +/- spin-orbit split columns).
The non-relativistic columns (TE_nonrel and the bare nlj columns) are
dropped.  Empty cells become null.
"""
import argparse
import csv
import json
import os

# Relativistic orbital columns, in the requested order.
ORBITALS = [
    "1s+",
    "2s+", "2p-", "2p+",
    "3s+", "3p-", "3p+", "3d-", "3d+",
    "4s+", "4p-", "4p+", "4d-", "4d+", "4f-", "4f+",
    "5s+", "5p-", "5p+", "5d-", "5d+", "5f-", "5f+",
    "6s+", "6p-", "6p+", "6d-", "6d+",
    "7s+",
]
FIELDS = ["Z", "Conf", "TE_rel"] + ORBITALS


def cell(value):
    """Strip a CSV cell; return None if empty."""
    value = (value or "").strip()
    return value if value else None


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", default=os.path.join(here, "DHF_GS_Energies.csv"))
    ap.add_argument("-o", "--output", default=os.path.join(here, "DHF_GS_Energies_rel.json"))
    args = ap.parse_args()

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        # DictReader keys retain trailing spaces from the header, so normalize.
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        missing = [c for c in FIELDS if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Input is missing expected columns: {missing}")

        records = []
        for row in reader:
            row = {k.strip(): v for k, v in row.items()}
            rec = {}
            for field in FIELDS:
                raw = cell(row.get(field))
                if field == "Z":
                    rec[field] = int(raw) if raw is not None else None
                elif field == "Conf":
                    rec[field] = raw
                else:
                    rec[field] = float(raw) if raw is not None else None
            records.append(rec)

    with open(args.output, "w") as f:
        json.dump(records, f, indent=2)
        f.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
