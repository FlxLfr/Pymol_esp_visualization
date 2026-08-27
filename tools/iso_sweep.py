#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iso_sweep.py
============

Backs the choice of ``--iso`` with measurements instead of a gut feeling.

The isovalue is the only parameter of the workflow that shifts the *measured*
numbers and not just the picture. This script shows by how much: it evaluates
the same cube files at several isovalues and puts the sigma hole, the surface
extrema and the belt side by side.

The evaluation uses the same functions ``render_esp.py`` uses
(``esp_statistics``, ``shell_points``, ``local_extrema``,
``sigma_hole_interpolated``). The table therefore holds exactly the numbers the
workflow itself would report; there is no second, diverging implementation that
could drift apart later.

Call::

    cd tools
    python iso_sweep.py --folder ../sandbox/brombenzol

Result: a table on the console and ``iso_sweep_<folder>.csv`` next to the cube
files. The table in section 4.1 of the background document comes from this
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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "scripts"))

import render_esp as R                                          # noqa: E402
from xyzToCube import esp_statistics                            # noqa: E402
from constants import HARTREE_TO_KCAL                           # noqa: E402


def measure(dens, esp, atoms, origin, voxel, iso):
    """The statistics for one isovalue.

    Everything large stays local and is released at the end: at 251^3 several
    hundred MB hang on every pass, which would otherwise sit around until the
    end of the loop and topple the process on a 3 GB machine.
    """
    vmin, vmax, npts = esp_statistics(dens, esp, iso=iso)
    pos, vals = R.shell_points(dens, esp, origin, voxel, iso=iso)
    loc = R.local_extrema(pos, vals, atoms)
    point_based = (loc.get("halogens") or [{}])[0].get("sigma_max")

    ray = R.sigma_hole_interpolated(dens, esp, origin, voxel, atoms, iso=iso)
    for entry in loc.get("halogens", []):
        if entry["index"] in ray:
            entry.update(ray[entry["index"]])
    R.promote_primary(loc)
    h = (loc.get("halogens") or [{}])[0]

    row = {
        "iso_au": iso,
        "shell_points": npts,
        "VS_min_au": vmin,
        "VS_min_on": loc.get("vmin_atom"),
        "VS_max_au": vmax,
        "VS_max_on": loc.get("vmax_atom"),
        "sigma_hole_au": h.get("sigma_max"),
        "sigma_hole_pointbased_au": point_based,
        "sigma_angle_deg": h.get("sigma_angle"),
        "belt_min_au": h.get("belt_min"),
        "belt_points": h.get("belt_points"),
        "halogen": h.get("label"),
    }
    del pos, vals, loc, ray
    gc.collect()
    return row


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(
        description="Plot the isovalue against the measured statistics.")
    p.add_argument("--folder", default=os.path.join(here, "..", "sandbox",
                                                    "brombenzol"),
                   help="molecule folder with td.cube and tp.cube "
                        "(default: ../sandbox/brombenzol)")
    p.add_argument("--isos", type=float, nargs="+",
                   default=[0.0005, 0.0008, 0.0010, 0.0015, 0.0020, 0.0040],
                   help="isovalues to test, in a.u. "
                        "(default: 0.0005 0.0008 0.001 0.0015 0.002 0.004)")
    p.add_argument("--out", default=None,
                   help="CSV output (default: iso_sweep_<folder>.csv next to "
                        "the cube files)")
    args = p.parse_args(argv)

    folder = os.path.abspath(args.folder)
    name = os.path.basename(os.path.normpath(folder))
    out = args.out or os.path.join(folder, f"iso_sweep_{name}.csv")

    print(f"reading {os.path.join(folder, 'td.cube')} ...", flush=True)
    dens, atoms, origin, voxel = R.read_cube(os.path.join(folder, "td.cube"))
    gc.collect()
    print(f"reading {os.path.join(folder, 'tp.cube')} ...", flush=True)
    esp, _, _, _ = R.read_cube(os.path.join(folder, "tp.cube"))
    gc.collect()
    if dens.shape != esp.shape:
        raise SystemExit("Density and ESP cube are on different grids.")
    print(f"grid {'x'.join(map(str, dens.shape))}, {len(atoms)} atoms\n",
          flush=True)

    rows = [measure(dens, esp, atoms, origin, voxel, iso) for iso in args.isos]

    hdr = (f"{'rho':>8} {'shell':>7} {'VS,min':>9} {'VS,max':>9} "
           f"{'sigma':>9} {'kcal':>7} {'belt':>9}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        sig, belt = r["sigma_hole_au"], r["belt_min_au"]
        print(f"{r['iso_au']:>8.4f} {r['shell_points']:>7} "
              f"{r['VS_min_au']:>+9.4f} {r['VS_max_au']:>+9.4f} "
              + (f"{sig:>+9.4f} {sig * HARTREE_TO_KCAL:>7.1f} "
                 if sig is not None else f"{'--':>9} {'--':>7} ")
              + (f"{belt:>+9.4f}" if belt is not None else f"{'--':>9}"))

    fields = ["iso_au", "shell_points", "VS_min_au", "VS_min_on", "VS_max_au",
              "VS_max_on", "sigma_hole_au", "sigma_hole_kcal",
              "sigma_hole_pointbased_au", "sigma_angle_deg", "belt_min_au",
              "belt_points", "halogen"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            r = dict(r)
            r["sigma_hole_kcal"] = ("" if r["sigma_hole_au"] is None else
                                    f"{r['sigma_hole_au'] * HARTREE_TO_KCAL:.2f}")
            r["iso_au"] = f"{r['iso_au']:.4f}"
            for k in ("VS_min_au", "VS_max_au", "sigma_hole_au",
                      "sigma_hole_pointbased_au", "belt_min_au"):
                r[k] = "" if r[k] is None else f"{r[k]:.5f}"
            r["sigma_angle_deg"] = ("" if r["sigma_angle_deg"] is None
                                    else f"{r['sigma_angle_deg']:.1f}")
            w.writerow(r)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
