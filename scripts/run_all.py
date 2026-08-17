#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py
==========

Batch driver for the ESP visualisation workflow.

Walks a directory tree, and for every molecule folder it finds:

  1. converts the Turbomole ``pointval`` grids to Gaussian cube (if the cube
     files are not there yet),
  2. renders the standard set of ESP images,
  3. records V_S,min and V_S,max.

Finally it writes ``summary.csv`` with the surface ESP statistics of every
molecule, and tells you which colour-scale range covers all of them.


A molecule folder is any directory that contains either

    td.xyz  +  tp.xyz          (Turbomole output, will be converted)
or  td.cube +  tp.cube         (already converted)

plus a structure file (``.mol``, ``.sdf``, ``.pdb`` or ``.xyz``).


Typical use
-----------
First pass, every molecule on its own automatic colour scale::

    python run_all.py --root ../examples --stride 2

Second pass, one common scale so the images are directly comparable::

    python run_all.py --root ../examples --stride 2 --esp-range 0.035

Or let the script do both in one go -- it runs the automatic pass, takes the
largest range it saw, and re-renders everything with it::

    python run_all.py --root ../examples --stride 2 --two-pass

``--two-pass`` is the recommended way to produce a comparable figure set:
you get the per-molecule maximum contrast images *and* the comparable ones,
without having to copy a number by hand.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import xyzToCube                                    # noqa: E402
import render_esp                                   # noqa: E402


HARTREE_TO_KCAL = 627.5095
HARTREE_TO_KJ = 2625.4996

STRUCT_EXT = (".mol", ".sdf", ".pdb", ".xyz")
GRID_NAMES = ("td", "tp")


# ----------------------------------------------------------------------------
# Discovery
# ----------------------------------------------------------------------------

def find_structure(folder, exclude):
    """First structure file in ``folder`` that is not one of the grid files."""
    excl = {os.path.abspath(p) for p in exclude if p}
    # .mol/.sdf/.pdb preferred over .xyz: they carry bond information, and a
    # bare .xyz can collide with the Turbomole grid files of the same suffix.
    for ext in STRUCT_EXT:
        for name in sorted(os.listdir(folder)):
            if not name.lower().endswith(ext):
                continue
            path = os.path.join(folder, name)
            if os.path.abspath(path) in excl:
                continue
            if os.path.splitext(name)[0] in GRID_NAMES:
                continue                      # td.xyz / tp.xyz are data, not structure
            return path
    return None


def discover(root):
    """All molecule folders below ``root``, sorted by name."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in ("images", "__pycache__", ".git")]
        names = set(filenames)
        has_raw = {"td.xyz", "tp.xyz"} <= names
        has_cube = {"td.cube", "tp.cube"} <= names
        if not (has_raw or has_cube):
            continue
        struct = find_structure(
            dirpath, exclude=[os.path.join(dirpath, "td.xyz"),
                              os.path.join(dirpath, "tp.xyz")])
        if struct is None:
            print(f"  ! {dirpath}: grids found but no structure file "
                  f"({'/'.join(STRUCT_EXT)}) - skipped")
            continue
        found.append({"dir": dirpath, "struct": struct,
                      "has_raw": has_raw, "has_cube": has_cube})
    return sorted(found, key=lambda e: e["dir"])


# ----------------------------------------------------------------------------
# Steps
# ----------------------------------------------------------------------------

def convert(entry, stride, struct_unit, force=False):
    """Turbomole grids -> cube, unless the cube files already exist."""
    folder = entry["dir"]
    out = {}
    for tag in GRID_NAMES:
        raw = os.path.join(folder, f"{tag}.xyz")
        cube = os.path.join(folder, f"{tag}.cube")
        out[tag] = cube
        if os.path.exists(cube) and not force:
            print(f"    {tag}.cube exists, skipping conversion")
            continue
        if not os.path.exists(raw):
            raise SystemExit(f"{folder}: neither {tag}.cube nor {tag}.xyz")
        print(f"    converting {tag}.xyz -> {tag}.cube (stride {stride})")
        atoms = xyzToCube.read_structure(entry["struct"], unit=struct_unit)
        info, data = xyzToCube.read_values(raw, verbose=False)
        xyzToCube.write_cube(cube, info, data, atoms, stride=stride,
                             comment=f"{info['quantity']} from {tag}.xyz")
    return out


def render(entry, cubes, esp_range, iso, transparency, backgrounds,
           width, height, dpi, buffer, prefix=None):
    folder = entry["dir"]
    args = types.SimpleNamespace(
        density=cubes["td"],
        esp=cubes["tp"],
        struct=entry["struct"],
        prefix=prefix or os.path.basename(os.path.normpath(folder)),
        outdir=os.path.join(folder, "images"),
        iso=iso,
        esp_range=esp_range,
        transparency=transparency,
        backgrounds=backgrounds,
        views=None,
        width=width,
        height=height,
        dpi=dpi,
        buffer=buffer,
    )
    return render_esp.render_all(args)


def write_summary(path, rows, common_range=None):
    fields = ["molecule", "structure", "grid", "iso_au",
              "shell_points", "VS_min_au", "VS_max_au",
              "VS_min_kcal", "VS_max_kcal", "VS_min_kJ", "VS_max_kJ",
              "esp_range_used_au", "esp_range_mode"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({
                "molecule": r["prefix"],
                "structure": os.path.basename(r["struct"]),
                "grid": r["grid"],
                "iso_au": r["iso"],
                "shell_points": r["shell_points"],
                "VS_min_au": f"{r['vmin']:.5f}",
                "VS_max_au": f"{r['vmax']:.5f}",
                "VS_min_kcal": f"{r['vmin'] * HARTREE_TO_KCAL:.2f}",
                "VS_max_kcal": f"{r['vmax'] * HARTREE_TO_KCAL:.2f}",
                "VS_min_kJ": f"{r['vmin'] * HARTREE_TO_KJ:.1f}",
                "VS_max_kJ": f"{r['vmax'] * HARTREE_TO_KJ:.1f}",
                "esp_range_used_au": f"{r['esp_range']:.4f}",
                "esp_range_mode": r["esp_range_mode"],
            })
        if common_range is not None:
            fh.write(f"# common colour scale covering all molecules: "
                     f"+/- {common_range:.4f} a.u.\n")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Batch conversion and rendering of molecular ESP data.")
    p.add_argument("--root", default="examples",
                   help="directory tree to search for molecule folders")
    p.add_argument("--stride", type=int, default=2,
                   help="grid decimation during conversion (default 2)")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom")
    p.add_argument("--force-convert", action="store_true",
                   help="rewrite cube files even if they already exist")
    p.add_argument("--esp-range", default="auto",
                   help="'auto' (per molecule) or a fixed value in a.u.")
    p.add_argument("--two-pass", action="store_true",
                   help="auto pass first, then re-render everything with the "
                        "largest range found")
    p.add_argument("--iso", type=float, default=0.001)
    p.add_argument("--transparency", type=float, default=0.15)
    p.add_argument("--backgrounds", nargs="+", default=["white"])
    p.add_argument("--width", type=int, default=2000)
    p.add_argument("--height", type=int, default=1600)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--buffer", type=float, default=2.4)
    p.add_argument("--summary", default=None,
                   help="path of the CSV summary (default <root>/summary.csv)")
    args = p.parse_args(argv)

    print("=" * 70)
    print("run_all.py - batch ESP visualisation")
    print("=" * 70)

    entries = discover(args.root)
    if not entries:
        raise SystemExit(f"No molecule folders found below '{args.root}'.")
    print(f"{len(entries)} molecule folder(s) found:")
    for e in entries:
        print(f"  - {e['dir']}  (structure: {os.path.basename(e['struct'])})")

    rows = []
    for e in entries:
        print(f"\n[{os.path.basename(os.path.normpath(e['dir']))}]")
        cubes = convert(e, args.stride, args.struct_unit, args.force_convert)
        res = render(e, cubes, args.esp_range, args.iso, args.transparency,
                     args.backgrounds, args.width, args.height, args.dpi,
                     args.buffer)
        rows.append(res)

    common = max(r["esp_range"] for r in rows)

    if args.two_pass and args.esp_range == "auto" and len(rows) > 1:
        print(f"\nSecond pass with a common colour scale of "
              f"+/- {common:.4f} a.u.")
        rows = []
        for e in entries:
            print(f"\n[{os.path.basename(os.path.normpath(e['dir']))}]")
            cubes = {t: os.path.join(e["dir"], f"{t}.cube")
                     for t in GRID_NAMES}
            res = render(e, cubes, common, args.iso, args.transparency,
                         args.backgrounds, args.width, args.height, args.dpi,
                         args.buffer)
            rows.append(res)
    elif args.two_pass and len(rows) <= 1:
        print("\n(--two-pass skipped: a single molecule needs no common scale)")

    summary = args.summary or os.path.join(args.root, "summary.csv")
    write_summary(summary, rows, common_range=common)

    print("\n" + "-" * 70)
    print(f"{'molecule':<24}{'V_S,min':>12}{'V_S,max':>12}"
          f"{'range used':>14}")
    print("-" * 70)
    for r in rows:
        print(f"{r['prefix']:<24}{r['vmin']:>+12.4f}{r['vmax']:>+12.4f}"
              f"{r['esp_range']:>14.4f}")
    print("-" * 70)
    print(f"Common scale covering all molecules: +/- {common:.4f} a.u.")
    print(f"Summary written to {summary}")
    return 0


if __name__ != "run_all":
    _argv = sys.argv[1:]
    if "--" in _argv:
        _argv = _argv[_argv.index("--") + 1:]
    else:
        for _i, _a in enumerate(_argv):
            if _a.endswith("run_all.py"):
                _argv = _argv[_i + 1:]
                break
    main(_argv)
