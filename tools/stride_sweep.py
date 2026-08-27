#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stride_sweep.py
===============

Backs the choice of ``--stride`` with measurements instead of a gut feeling.

For every stride it does the same thing ``xyzToCube.py --stride N`` does to the
cube files - keep every N-th grid point per axis, stretch the voxel vectors
accordingly - and then determines the statistics with the same functions
``render_esp.py`` uses. The table therefore holds exactly the numbers the
workflow itself would report; there is no second, diverging implementation.

Why not reconvert the pointval files N times: the result would be the same
(write_cube decimates in exactly this way), but every pass would have to parse
1.25 GB of text. From the finished cubes the same test is done in seconds.

Call::

    cd tools
    python stride_sweep.py --folder ../sandbox/brombenzol

Result: a table on the console and ``stride_sweep_<folder>.csv`` next to the
cube files. The table in section 4.2 of the background document comes from this
script.

Only numpy and the scripts in ../scripts are needed - no PyMOL: render_esp.py
imports that only inside ensure_pymol().
"""

from __future__ import annotations

import argparse
import csv
import gc
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "scripts"))

import render_esp as R                                          # noqa: E402
from xyzToCube import esp_statistics                            # noqa: E402
from constants import HARTREE_TO_KCAL                           # noqa: E402


def cube_bytes(n, natoms, newline=2):
    """Size of the cube file write_cube would produce for this grid.

    Computed rather than measured, because the file for stride 1 alone would be
    about 200 MB of write load - not appropriate for one table column.

    write_cube writes the nz values per (x, y) pair with ``%13.5E`` and breaks
    the line after every six; the header is 2 comment lines, 1 origin line,
    3 axis lines and one line per atom. ``newline`` is the length of the line
    break the text mode turns that into: 2 on Windows (CRLF), 1 on
    Linux/macOS.

    The data part is exact. In the header only the two comment lines vary,
    because the name of the source file appears in them; they are estimated at
    35 characters. Checked against sandbox/brombenzol/td.cube (251^3, 12 atoms,
    written on Windows): computed 210_865_313, measured 210_865_313 bytes. With
    other file names the header differs by a few bytes - immaterial at 201 MB.
    """
    nx, ny, nz = n
    lines_per_row = -(-nz // 6)                       # rounded up
    body = nx * ny * (nz * 13 + lines_per_row * newline)
    # Header: 2 comment lines, 1 atom-count line, 3 axis lines, 1 line per
    # atom. Line lengths from the format strings in write_cube: origin and axis
    # lines 44 characters, atom lines 57. The two comment lines depend on the
    # file name and are estimated at 35 characters.
    header = (4 * (44 + newline) + natoms * (57 + newline)
              + 2 * (35 + newline))
    return body + header


def measure(dens, esp, atoms, origin, voxel, stride, iso):
    """The statistics for a grid decimated by ``stride``."""
    d = np.ascontiguousarray(dens[::stride, ::stride, ::stride])
    e = np.ascontiguousarray(esp[::stride, ::stride, ::stride])
    v = voxel * stride

    vmin, vmax, npts = esp_statistics(d, e, iso=iso)
    pos, vals = R.shell_points(d, e, origin, v, iso=iso)
    loc = R.local_extrema(pos, vals, atoms)
    point_based = (loc.get("halogens") or [{}])[0].get("sigma_max")

    ray = R.sigma_hole_interpolated(d, e, origin, v, atoms, iso=iso)
    for entry in loc.get("halogens", []):
        if entry["index"] in ray:
            entry.update(ray[entry["index"]])
    R.promote_primary(loc)
    h = (loc.get("halogens") or [{}])[0]

    row = {
        "stride": stride,
        "grid": "x".join(str(s) for s in d.shape),
        "spacing_bohr": float(np.max(np.abs(np.diag(v)))),
        "cube_bytes": cube_bytes(d.shape, len(atoms)),
        "shell_points": npts,
        "VS_min_au": vmin,
        "VS_min_on": loc.get("vmin_atom"),
        "VS_max_au": vmax,
        "VS_max_on": loc.get("vmax_atom"),
        "sigma_hole_au": h.get("sigma_max"),
        "sigma_hole_pointbased_au": point_based,
        "belt_min_au": h.get("belt_min"),
        "belt_points": h.get("belt_points"),
        "halogen": h.get("label"),
    }
    del d, e, pos, vals, loc, ray
    gc.collect()
    return row


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(
        description="Plot the grid resolution (--stride) against the measured "
                    "statistics.")
    p.add_argument("--folder", default=os.path.join(here, "..", "sandbox",
                                                    "brombenzol"),
                   help="molecule folder with td.cube and tp.cube "
                        "(default: ../sandbox/brombenzol)")
    p.add_argument("--strides", type=int, nargs="+", default=[1, 2, 3, 4, 6, 8],
                   help="stride values to test (default: 1 2 3 4 6 8)")
    p.add_argument("--iso", type=float, default=0.001,
                   help="isovalue, fixed across all strides (default: 0.001)")
    p.add_argument("--out", default=None,
                   help="CSV output (default: stride_sweep_<folder>.csv "
                        "neben den Cube-Dateien)")
    args = p.parse_args(argv)

    folder = os.path.abspath(args.folder)
    name = os.path.basename(os.path.normpath(folder))
    out = args.out or os.path.join(folder, f"stride_sweep_{name}.csv")

    print(f"reading {os.path.join(folder, 'td.cube')} ...", flush=True)
    dens, atoms, origin, voxel = R.read_cube(os.path.join(folder, "td.cube"))
    gc.collect()
    print(f"reading {os.path.join(folder, 'tp.cube')} ...", flush=True)
    esp, _, _, _ = R.read_cube(os.path.join(folder, "tp.cube"))
    gc.collect()
    if dens.shape != esp.shape:
        raise SystemExit("Density and ESP cube are on different grids.")
    print(f"grid {'x'.join(map(str, dens.shape))}, {len(atoms)} atoms, "
          f"isovalue {args.iso}\n", flush=True)

    rows = [measure(dens, esp, atoms, origin, voxel, s, args.iso)
            for s in args.strides]

    hdr = (f"{'stride':>6} {'grid':>12} {'d/Bohr':>7} {'cubes':>9} "
           f"{'shell':>7} {'VS,min':>9} {'VS,max':>9} {'sigma':>9} "
           f"{'kcal':>7} {'belt':>9}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        sig = r["sigma_hole_au"]
        belt = r["belt_min_au"]
        print(f"{r['stride']:>6} {r['grid']:>12} {r['spacing_bohr']:>7.2f} "
              f"{r['cube_bytes'] / 1024**2:>7.1f} MB {r['shell_points']:>7} "
              f"{r['VS_min_au']:>+9.4f} {r['VS_max_au']:>+9.4f} "
              f"{sig:>+9.4f} {sig * HARTREE_TO_KCAL:>7.1f} "
              + (f"{belt:>+9.4f}" if belt is not None else f"{'--':>9}"))

    fields = ["stride", "grid", "spacing_bohr", "cube_bytes", "shell_points",
              "VS_min_au", "VS_min_on", "VS_max_au", "VS_max_on",
              "sigma_hole_au", "sigma_hole_kcal", "sigma_hole_pointbased_au",
              "belt_min_au", "belt_points", "halogen"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            r = dict(r)
            r["spacing_bohr"] = f"{r['spacing_bohr']:.4f}"
            for k in ("VS_min_au", "VS_max_au", "sigma_hole_au",
                      "sigma_hole_pointbased_au", "belt_min_au"):
                r[k] = "" if r[k] is None else f"{r[k]:.5f}"
            r["sigma_hole_kcal"] = ("" if not r["sigma_hole_au"] else
                                    f"{float(r['sigma_hole_au']) * HARTREE_TO_KCAL:.2f}")
            w.writerow(r)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
