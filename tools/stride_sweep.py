#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stride_sweep.py
===============

Belegt die Wahl von ``--stride`` mit Messwerten statt mit einem Bauchgefuehl.

Fuer jeden Stride wird dasselbe gemacht, was ``xyzToCube.py --stride N`` mit den
Cube-Dateien tut - jeden N-ten Gitterpunkt pro Achse behalten, Voxelvektoren
entsprechend strecken - und anschliessend werden die Kennwerte mit denselben
Funktionen bestimmt, die auch ``render_esp.py`` benutzt. Dadurch stehen in der
Tabelle exakt die Zahlen, die der Workflow selbst ausgeben wuerde; es gibt keine
zweite, abweichende Implementierung.

Warum nicht die pointval-Dateien N-mal neu konvertieren: das Ergebnis waere
dasselbe (write_cube dezimiert genau so), aber jeder Durchlauf muesste 1.25 GB
Text parsen. Aus den fertigen Cubes ist derselbe Test in Sekunden erledigt.

Aufruf::

    cd tools
    python stride_sweep.py --folder ../sandbox/brombenzol

Ergebnis: Tabelle auf der Konsole und ``stride_sweep_<ordner>.csv`` daneben.
Die Tabelle in Abschnitt 4.2 der Hintergrunddokumentation stammt aus diesem
Skript.

Gebraucht werden nur numpy und die Skripte in ../scripts - kein PyMOL: das wird
in render_esp.py erst in ensure_pymol() importiert.
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
    """Groesse der Cube-Datei, die write_cube fuer dieses Gitter schreiben wuerde.

    Rechnung statt Messung, weil allein die Datei fuer Stride 1 rund 200 MB
    Schreiblast waere - fuer eine Tabellenspalte ist das nicht angemessen.

    write_cube schreibt pro (x, y)-Paar die nz Werte mit ``%13.5E`` und bricht
    nach je sechs Werten um; Kopf sind 2 Kommentarzeilen, 1 Ursprungszeile,
    3 Achsenzeilen und je Atom eine Zeile. ``newline`` ist die Laenge des
    Zeilenumbruchs, den der Textmodus daraus macht: 2 unter Windows (CRLF),
    1 unter Linux/macOS.

    Der Datenteil ist exakt. Im Kopf sind nur die zwei Kommentarzeilen
    variabel, weil der Name der Quelldatei darin steht; sie werden mit
    35 Zeichen veranschlagt. Geprueft an sandbox/brombenzol/td.cube
    (251^3, 12 Atome, unter Windows geschrieben): berechnet 210_865_313,
    gemessen 210_865_313 Byte. Bei anderen Dateinamen weicht der Kopf um
    einige Byte ab - auf 201 MB ohne Bedeutung.
    """
    nx, ny, nz = n
    lines_per_row = -(-nz // 6)                       # aufgerundet
    body = nx * ny * (nz * 13 + lines_per_row * newline)
    # Kopf: 2 Kommentarzeilen, 1 Atomzahlzeile, 3 Achsenzeilen, je Atom 1 Zeile
    # Zeilenlaengen aus den Formatstrings in write_cube: Ursprungs- und
    # Achsenzeilen 44 Zeichen, Atomzeilen 57. Die zwei Kommentarzeilen haengen
    # vom Dateinamen ab und werden mit 35 Zeichen veranschlagt.
    header = (4 * (44 + newline) + natoms * (57 + newline)
              + 2 * (35 + newline))
    return body + header


def measure(dens, esp, atoms, origin, voxel, stride, iso):
    """Kennwerte fuer ein um ``stride`` dezimiertes Gitter."""
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
        description="Gitteraufloesung (--stride) gegen die gemessenen "
                    "Kennwerte auftragen.")
    p.add_argument("--folder", default=os.path.join(here, "..", "sandbox",
                                                    "brombenzol"),
                   help="Molekuelordner mit td.cube und tp.cube "
                        "(Standard: ../sandbox/brombenzol)")
    p.add_argument("--strides", type=int, nargs="+", default=[1, 2, 3, 4, 6, 8],
                   help="zu testende Stride-Werte (Standard: 1 2 3 4 6 8)")
    p.add_argument("--iso", type=float, default=0.001,
                   help="Isowert, ueber alle Strides fest (Standard: 0.001)")
    p.add_argument("--out", default=None,
                   help="CSV-Ausgabe (Standard: stride_sweep_<ordner>.csv "
                        "neben den Cube-Dateien)")
    args = p.parse_args(argv)

    folder = os.path.abspath(args.folder)
    name = os.path.basename(os.path.normpath(folder))
    out = args.out or os.path.join(folder, f"stride_sweep_{name}.csv")

    print(f"lese {os.path.join(folder, 'td.cube')} ...", flush=True)
    dens, atoms, origin, voxel = R.read_cube(os.path.join(folder, "td.cube"))
    gc.collect()
    print(f"lese {os.path.join(folder, 'tp.cube')} ...", flush=True)
    esp, _, _, _ = R.read_cube(os.path.join(folder, "tp.cube"))
    gc.collect()
    if dens.shape != esp.shape:
        raise SystemExit("Dichte- und ESP-Cube haben unterschiedliche Gitter.")
    print(f"Gitter {'x'.join(map(str, dens.shape))}, {len(atoms)} Atome, "
          f"Isowert {args.iso}\n", flush=True)

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
