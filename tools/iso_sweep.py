#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iso_sweep.py
============

Belegt die Wahl von ``--iso`` mit Messwerten statt mit einem Bauchgefuehl.

Der Isowert ist der einzige Parameter des Workflows, der die *gemessenen* Zahlen
verschiebt und nicht nur das Bild. Dieses Skript zeigt, wie stark: es wertet
dieselben Cube-Dateien bei mehreren Isowerten aus und stellt sigma-Loch,
Oberflaechenextrema und Guertel nebeneinander.

Ausgewertet wird mit denselben Funktionen, die auch ``render_esp.py`` benutzt
(``esp_statistics``, ``shell_points``, ``local_extrema``,
``sigma_hole_interpolated``). Dadurch stehen in der Tabelle exakt die Zahlen,
die der Workflow selbst ausgeben wuerde; es gibt keine zweite, abweichende
Implementierung, die spaeter auseinanderlaufen koennte.

Aufruf::

    cd tools
    python iso_sweep.py --folder ../sandbox/brombenzol

Ergebnis: Tabelle auf der Konsole und ``iso_sweep_<ordner>.csv`` daneben.
Die Tabelle in Abschnitt 4.1 der Hintergrunddokumentation stammt aus diesem
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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "scripts"))

import render_esp as R                                          # noqa: E402
from xyzToCube import esp_statistics                            # noqa: E402
from constants import HARTREE_TO_KCAL                           # noqa: E402


def measure(dens, esp, atoms, origin, voxel, iso):
    """Kennwerte fuer einen Isowert.

    Alles Grosse bleibt lokal und wird am Ende freigegeben: bei 251^3 haengen an
    jedem Durchlauf mehrere hundert MB, die sonst bis zum Schleifenende liegen
    bleiben und den Prozess auf einer 3-GB-Maschine kippen.
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
        description="Isowert gegen die gemessenen Kennwerte auftragen.")
    p.add_argument("--folder", default=os.path.join(here, "..", "sandbox",
                                                    "brombenzol"),
                   help="Molekuelordner mit td.cube und tp.cube "
                        "(Standard: ../sandbox/brombenzol)")
    p.add_argument("--isos", type=float, nargs="+",
                   default=[0.0005, 0.0008, 0.0010, 0.0015, 0.0020, 0.0040],
                   help="zu testende Isowerte in a.u. "
                        "(Standard: 0.0005 0.0008 0.001 0.0015 0.002 0.004)")
    p.add_argument("--out", default=None,
                   help="CSV-Ausgabe (Standard: iso_sweep_<ordner>.csv "
                        "neben den Cube-Dateien)")
    args = p.parse_args(argv)

    folder = os.path.abspath(args.folder)
    name = os.path.basename(os.path.normpath(folder))
    out = args.out or os.path.join(folder, f"iso_sweep_{name}.csv")

    print(f"lese {os.path.join(folder, 'td.cube')} ...", flush=True)
    dens, atoms, origin, voxel = R.read_cube(os.path.join(folder, "td.cube"))
    gc.collect()
    print(f"lese {os.path.join(folder, 'tp.cube')} ...", flush=True)
    esp, _, _, _ = R.read_cube(os.path.join(folder, "tp.cube"))
    gc.collect()
    if dens.shape != esp.shape:
        raise SystemExit("Dichte- und ESP-Cube haben unterschiedliche Gitter.")
    print(f"Gitter {'x'.join(map(str, dens.shape))}, {len(atoms)} Atome\n",
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
