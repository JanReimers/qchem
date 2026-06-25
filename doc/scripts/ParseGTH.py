#!/usr/bin/env python3
"""
Parse the CP2K GTH/HGH pseudopotential database (BasisSetData/GTH_POTENTIALS.txt)
into hierarchical JSON for the C++ reader.

Source: CP2K data/GTH_POTENTIALS (Goedecker-Teter-Hutter / Hartwigsen-Goedecker-Hutter
        analytic norm-conserving pseudopotentials).
        Goedecker, Teter, Hutter, PRB 54, 1703 (1996); Hartwigsen, Goedecker, Hutter,
        PRB 58, 3641 (1998).

File format (one block per element/functional/valence; '#' comments, blank lines ignored):

    Symbol  Name  [Aliases...]            # e.g. "Si GTH-PADE-q4 GTH-LDA-q4 GTH-PADE GTH-LDA"
    n_elec(s) n_elec(p) n_elec(d) ...     # valence config; sum = q = Zion
    r_loc  nexp  C1 ... C_nexp            # local part (nexp may be 0)
    nprj                                  # number of non-local channels (l = 0,1,2,...; may be 0)
    r_l  nprj_l  <upper triangle of the nprj_l x nprj_l h-matrix, row-major,
                  spilling onto continuation lines>
    ...                                   # one such group per channel

The h-matrix upper triangle holds nprj_l*(nprj_l+1)/2 numbers; we read them as a flat
token stream (continuation lines have no leading r_l/nprj_l) and mirror to a full
symmetric matrix.

Output JSON, keyed for O(1) lookup by (element, functional, q):

    { "Si": { "LDA": { "default": "4",
                       "4": { "z_ion": 4, "n_elec": [2,2], "r_loc": 0.44,
                              "c": [-7.336...],
                              "channels": [ { "l": 0, "r": 0.4227,
                                              "h": [[5.9069,-1.2619],[-1.2619,3.2582]] },
                                            { "l": 1, "r": 0.4843, "h": [[2.7270]] } ] } } } }

Usage:  python3 doc/scripts/ParseGTH.py \
            BasisSetData/GTH_POTENTIALS.txt \
            src/BasisSet/Lattice_3D/Data/gth_potentials.json
"""

import json
import re
import sys
from pathlib import Path

# CP2K calls the LDA (Perdew-Zunger / Pade) functional both "PADE" and "LDA"; canonicalise to one key.
FUNCTIONAL_ALIASES = {"PADE": "LDA"}


def functional_of(name):
    """('GTH-PADE-q4') -> ('LDA', '4');  ('GTH-PBE') -> ('PBE', None)."""
    m = re.match(r"GTH-([A-Za-z0-9]+?)(?:-q(\d+))?$", name)
    if not m:
        return None, None
    func, q = m.group(1).upper(), m.group(2)
    return FUNCTIONAL_ALIASES.get(func, func), q


def read_blocks(lines):
    """Split the file into (header_tokens, [data_lines]); a header starts with an element symbol."""
    blocks = []
    header, data = None, []
    for line in lines:
        s = line.split("#", 1)[0].rstrip()
        if not s.strip():
            continue
        if re.match(r"^[A-Z][a-z]?\s+\S", s):          # element symbol begins a new block
            if header is not None:
                blocks.append((header, data))
            header, data = s.split(), []
        else:
            data.append(s)
    if header is not None:
        blocks.append((header, data))
    return blocks


def parse_block(header, data):
    """Return (symbol, names, record-dict) for one pseudopotential block."""
    symbol, names = header[0], header[1:]

    n_elec = [int(x) for x in data[0].split()]         # line 2: valence config
    z_ion = sum(n_elec)

    tok = " ".join(data[1:]).split()                   # flat number stream for the rest
    i = 0

    def take_float():
        nonlocal i
        v = float(tok[i]); i += 1; return v

    def take_int():
        nonlocal i
        v = int(tok[i]); i += 1; return v

    r_loc = take_float()
    nexp = take_int()
    c = [take_float() for _ in range(nexp)]

    channels = []
    nprj = take_int()
    for l in range(nprj):
        r_l = take_float()
        np_l = take_int()
        h = [[0.0] * np_l for _ in range(np_l)]        # mirror the upper triangle to full symmetric
        for a in range(np_l):
            for b in range(a, np_l):
                v = take_float()
                h[a][b] = h[b][a] = v
        channels.append({"l": l, "r": r_l, "h": h})

    assert i == len(tok), f"{symbol}: {len(tok)-i} leftover tokens (format mismatch)"
    return symbol, names, {"z_ion": z_ion, "n_elec": n_elec, "r_loc": r_loc, "c": c, "channels": channels}


def build_table(blocks):
    table = {}
    for header, data in blocks:
        symbol, names, rec = parse_block(header, data)
        q = str(rec["z_ion"])
        for name in names:
            func, qn = functional_of(name)
            if func is None:
                continue
            slot = table.setdefault(symbol, {}).setdefault(func, {})
            if qn is None:                              # the bare (no -q) alias marks the default valence
                slot["default"] = q
            else:
                assert qn == q, f"{symbol} {name}: name q{qn} != sum(n_elec)={q}"
                slot[q] = rec
    return table


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("BasisSetData/GTH_POTENTIALS.txt")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("src/BasisSet/Lattice_3D/Data/gth_potentials.json")

    blocks = read_blocks(src.read_text().splitlines())
    table = build_table(blocks)

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(table, indent=1, sort_keys=True))

    n_elem = len(table)
    n_pp = sum(len(q) - ("default" in fns) for el in table.values() for fns in el.values() for q in [fns])
    print(f"Parsed {len(blocks)} blocks -> {n_elem} elements into {dst}")
    if "Si" in table and "LDA" in table["Si"]:
        print("Si/LDA/q4:", json.dumps(table["Si"]["LDA"].get("4"), indent=1))


if __name__ == "__main__":
    main()
