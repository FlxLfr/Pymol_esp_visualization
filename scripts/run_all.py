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

plus a structure file (``.mol``, ``.sdf`` or ``.xyz``).


Typical use
-----------
First pass, every molecule on its own automatic colour scale::

    python run_all.py --root ../sandbox

Second pass, one common scale so the images are directly comparable::

    python run_all.py --root ../sandbox --esp-range 0.035

Or let the script do both in one go -- it runs the automatic pass, takes the
largest range it saw, and re-renders everything with it::

    python run_all.py --root ../sandbox --two-pass

Conversion runs at full grid resolution by default, like xyzToCube.py. For a
quicker first look, --stride 2 keeps every second point per axis: eight times
smaller cubes, visually identical images, and still enough grid points in the
sigma-hole cap. Note that --stride only takes effect while cube files are being
*created* -- if td.cube and tp.cube already exist they are reused as they are.

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
import ansi                                         # noqa: E402
from constants import HARTREE_TO_KCAL, HARTREE_TO_KJ  # noqa: E402

STRUCT_EXT = (".mol", ".sdf", ".xyz")
GRID_NAMES = ("td", "tp")

# reference/ liegt neben scripts/, nicht im aktuellen Arbeitsverzeichnis
DEFAULT_ROOT = os.path.normpath(os.path.join(_HERE, "..", "reference"))


# ----------------------------------------------------------------------------
# Discovery
# ----------------------------------------------------------------------------

def find_structure(folder, exclude):
    """First structure file in ``folder`` that is not one of the grid files."""
    excl = {os.path.abspath(p) for p in exclude if p}
    # .mol/.sdf preferred over .xyz: they carry bond information, and a
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


def find_cubes(folder, filenames):
    """Existing cube files in ``folder``, or None if not both are present.

    Accepts the decimated ``td_demo.cube`` / ``tp_demo.cube`` of the reference
    dataset as well, so that ``python run_all.py`` works as a smoke test on a
    fresh clone.
    """
    cubes = {}
    for tag in GRID_NAMES:
        for candidate in (f"{tag}.cube", f"{tag}_demo.cube"):
            if candidate in filenames:
                cubes[tag] = os.path.join(folder, candidate)
                break
    return cubes if len(cubes) == len(GRID_NAMES) else None


def discover(root):
    """All molecule folders below ``root``, sorted by name."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in ("images", "__pycache__", ".git")]
        names = set(filenames)
        has_raw = {"td.xyz", "tp.xyz"} <= names
        cubes = find_cubes(dirpath, names)
        if not (has_raw or cubes):
            continue
        struct = find_structure(
            dirpath, exclude=[os.path.join(dirpath, "td.xyz"),
                              os.path.join(dirpath, "tp.xyz")])
        if struct is None:
            print(f"  ! {dirpath}: grids found but no structure file "
                  f"({'/'.join(STRUCT_EXT)}) - skipped")
            continue
        found.append({"dir": dirpath, "struct": struct,
                      "has_raw": has_raw, "cubes": cubes})
    return sorted(found, key=lambda e: e["dir"])


# ----------------------------------------------------------------------------
# Steps
# ----------------------------------------------------------------------------

def convert(entry, stride, struct_unit, force=False):
    """Turbomole grids -> cube, unless the cube files already exist."""
    folder = entry["dir"]
    existing = entry.get("cubes") or {}
    out = {}
    for tag in GRID_NAMES:
        raw = os.path.join(folder, f"{tag}.xyz")
        cube = existing.get(tag) or os.path.join(folder, f"{tag}.cube")
        out[tag] = cube
        if tag in existing and not force:
            print(f"    {os.path.basename(cube)} exists, skipping conversion")
            continue
        if not os.path.exists(raw):
            raise SystemExit(f"{folder}: neither {tag}.cube nor {tag}.xyz")
        cube = os.path.join(folder, f"{tag}.cube")   # Neukonvertierung immer ohne _demo
        out[tag] = cube
        print(f"    converting {tag}.xyz -> {tag}.cube (stride {stride})")
        atoms = xyzToCube.read_structure(entry["struct"], unit=struct_unit)
        info, data = xyzToCube.read_values(raw, verbose=False)
        xyzToCube.write_cube(cube, info, data, atoms, stride=stride,
                             comment=f"{info['quantity']} from {tag}.xyz")
    return out


def write_scene(entry, cubes, esp_range, iso, transparency,
                filename="esp.pml", rainbow=False):
    """Schreibt ein PyMOL-Skript neben die Cube-Dateien.

    Damit hat jedes Molekuel nicht nur die fertigen Bilder, sondern auch eine
    interaktive Szene zum Drehen und Nachschauen - mit exakt der Farbskala, die
    auch fuer die Bilder verwendet wurde. Beim Zwei-Pass-Lauf wird die Datei im
    zweiten Durchgang mit der gemeinsamen Skala ueberschrieben, sie passt also
    immer zum zuletzt erzeugten Bildersatz.
    """
    folder = entry["dir"]
    path = os.path.join(folder, filename)
    xyzToCube.write_pymol_script(
        path,
        struct=os.path.basename(entry["struct"]),
        density_cube=os.path.basename(cubes["td"]),
        esp_cube=os.path.basename(cubes["tp"]),
        vmin=-esp_range, vmax=esp_range,
        iso=iso, transparency=transparency, rainbow=rainbow,
    )
    return path


def render(entry, cubes, esp_range, iso, transparency, backgrounds,
           width, height, dpi, buffer, prefix=None, images_dir="images",
           rainbow=False):
    folder = entry["dir"]
    args = types.SimpleNamespace(
        density=cubes["td"],
        esp=cubes["tp"],
        struct=entry["struct"],
        prefix=prefix or os.path.basename(os.path.normpath(folder)),
        outdir=os.path.join(folder, images_dir),
        iso=iso,
        esp_range=esp_range,
        transparency=transparency,
        backgrounds=backgrounds,
        views=None,
        width=width,
        height=height,
        dpi=dpi,
        buffer=buffer,
        rainbow=rainbow,
    )
    return render_esp.render_all(args)


def write_summary(path, rows, common_range=None):
    fields = ["molecule", "structure", "grid", "iso_au",
              "shell_points", "VS_min_au", "VS_max_au",
              "VS_min_kcal", "VS_max_kcal", "VS_min_kJ", "VS_max_kJ",
              "VS_max_on", "VS_min_on",
              # sigma_hole_* describe the STRONGEST sigma-hole; sigma_hole_on
              # says which atom that is, and sigma_holes_all lists every
              # halogen of the molecule (label:value in a.u., ';'-separated)
              # so nothing is lost for multi-halogen compounds.
              "halogen", "sigma_hole_on", "sigma_hole_au", "sigma_hole_kcal",
              "sigma_method", "belt_min_au", "belt_min_kcal",
              "n_halogens", "sigma_holes_all",
              "esp_range_used_au", "esp_range_mode", "colormap"]
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
                "VS_max_on": r.get("vmax_atom") or "",
                "VS_min_on": r.get("vmin_atom") or "",
                "halogen": r.get("halogen") or "",
                "sigma_hole_on": r.get("halogen_atom") or "",
                "n_halogens": r.get("n_halogens", 0),
                "sigma_holes_all": ";".join(
                    f"{e['label']}:{e['sigma_max']:.5f}"
                    for e in r.get("halogens", [])
                    if e.get("sigma_max") is not None),
                "sigma_hole_au": ("" if r.get("sigma_max") is None
                                  else f"{r['sigma_max']:.5f}"),
                "sigma_hole_kcal": ("" if r.get("sigma_max") is None
                                    else f"{r['sigma_max']*HARTREE_TO_KCAL:.2f}"),
                "sigma_method": r.get("sigma_method") or "",
                "belt_min_au": ("" if r.get("belt_min") is None
                                else f"{r['belt_min']:.5f}"),
                "belt_min_kcal": ("" if r.get("belt_min") is None
                                  else f"{r['belt_min']*HARTREE_TO_KCAL:.2f}"),
                "esp_range_used_au": f"{r['esp_range']:.4f}",
                "esp_range_mode": r["esp_range_mode"],
                "colormap": r.get("colormap", "redblue"),
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
    # Ohne --root wird der reference/-Ordner des Repositoriums benutzt, egal aus
    # welchem Verzeichnis das Skript aufgerufen wird. Das ist der Selbsttest.
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="directory tree to search for molecule folders "
                        "(default: the repository's reference/ folder)")
    p.add_argument("--only", nargs="+", metavar="NAME",
                   help="restrict the run to these molecule folders; simple "
                        "wildcards are allowed, e.g. --only paracetamol "
                        "'*benzol'")
    p.add_argument("--stride", type=int, default=1,
                   help="grid decimation during conversion: keep every n-th "
                        "point per axis (default 1 = full resolution, same as "
                        "xyzToCube.py). --stride 2 gives 8x smaller cubes and "
                        "is plenty for images; do not go coarser if you need "
                        "the sigma-hole value. Ignored when the cube files "
                        "already exist.")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom")
    p.add_argument("--force-convert", action="store_true",
                   help="rewrite cube files even if they already exist")
    p.add_argument("--esp-range", default="auto",
                   help="'auto' (per molecule) or a fixed value in a.u.")
    p.add_argument("--rainbow", action="store_true",
                   help="rainbow colour ramp instead of red-white-blue; "
                        "writes a separate <molecule>_rainbow_* image set "
                        "and esp_rainbow.pml, so the standard set is kept")
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
    p.add_argument("--images-dir", default=None,
                   help="name of the output folder inside each molecule "
                        "folder (default: 'images', or 'images_check' when "
                        "running the built-in reference set)")
    p.add_argument("--no-color", action="store_true",
                   help="plain output without ANSI colours (same effect as "
                        "setting the NO_COLOR environment variable)")
    p.add_argument("--summary", default=None,
                   help="path of the CSV summary (default <root>/summary.csv)")
    args = p.parse_args(argv)

    if args.no_color:
        ansi.disable()

    is_reference = os.path.abspath(args.root) == os.path.abspath(DEFAULT_ROOT)
    # Auch die PyMOL-Szene des Selbsttests darf nichts Committetes ueberschreiben.
    # Eigener Name, sonst ueberschreibt ein Regenbogenlauf die Szene des
    # rot-weiss-blauen Laufs - die Bilder liegen aus demselben Grund getrennt.
    _stem = "esp_check" if is_reference else "esp"
    pml_name = f"{_stem}{'_rainbow' if args.rainbow else ''}.pml"
    if args.images_dir is None:
        # Der Selbsttest darf die committeten Referenzbilder nicht ueberschreiben.
        args.images_dir = "images_check" if is_reference else "images"

    print("=" * 70)
    print("run_all.py - batch ESP visualisation")
    print("=" * 70)
    if is_reference:
        print("  Reference dataset (smoke test).")
        print(f"  Writing to '{args.images_dir}/' so that the committed "
              f"reference images stay untouched.")
        print("  Compare your output with reference/*/images/ - the colour "
              "scale and V_S values must match.")

    entries = discover(args.root)
    if args.only:
        import fnmatch
        keep = []
        for e in entries:
            base = os.path.basename(os.path.normpath(e["dir"]))
            if any(fnmatch.fnmatch(base.lower(), pat.lower())
                   for pat in args.only):
                keep.append(e)
        skipped = len(entries) - len(keep)
        entries = keep
        if skipped:
            print(f"  --only: {skipped} folder(s) skipped")
    if not entries:
        raise SystemExit(f"No molecule folders found below '{args.root}'"
                         + (" matching --only." if args.only else "."))
    print(f"{len(entries)} molecule folder(s) found:")
    for e in entries:
        print(f"  - {e['dir']}  (structure: {os.path.basename(e['struct'])})")

    rows = []
    for e in entries:
        name = os.path.basename(os.path.normpath(e["dir"]))
        print("\n" + ansi.paint(f"[{name}]", ansi.GREEN + ansi.BOLD))
        cubes = convert(e, args.stride, args.struct_unit, args.force_convert)
        res = render(e, cubes, args.esp_range, args.iso, args.transparency,
                     args.backgrounds, args.width, args.height, args.dpi,
                     args.buffer, images_dir=args.images_dir,
                     rainbow=args.rainbow)
        pml = write_scene(e, cubes, res["esp_range"], args.iso,
                          args.transparency, filename=pml_name,
                          rainbow=args.rainbow)
        print(f"    -> {pml}")
        rows.append(res)

    common = max(r["esp_range"] for r in rows)

    if args.two_pass and args.esp_range == "auto" and len(rows) > 1:
        print(f"\nSecond pass with a common colour scale of "
              f"+/- {common:.4f} a.u.")
        rows = []
        for e in entries:
            name = os.path.basename(os.path.normpath(e["dir"]))
            print("\n" + ansi.paint(f"[{name}]", ansi.GREEN + ansi.BOLD))
            cubes = convert(e, args.stride, args.struct_unit, force=False)
            res = render(e, cubes, common, args.iso, args.transparency,
                         args.backgrounds, args.width, args.height, args.dpi,
                         args.buffer, images_dir=args.images_dir,
                         rainbow=args.rainbow)
            pml = write_scene(e, cubes, common, args.iso, args.transparency,
                              filename=pml_name, rainbow=args.rainbow)
            print(f"    -> {pml}")
            rows.append(res)
    elif args.two_pass and len(rows) <= 1:
        print("\n(--two-pass skipped: a single molecule needs no common scale)")

    summary = args.summary or os.path.join(
        args.root, "summary_check.csv" if is_reference else "summary.csv")
    write_summary(summary, rows, common_range=common)

    print("\n" + "-" * 70)
    print(f"{'molecule':<20}{'V_S,min':>10}{'V_S,max':>10}{'(on)':>7}"
          f"{'sigma-hole':>12}{'range':>9}")
    print("-" * 70)
    for r in rows:
        # Platzhalter in derselben Spaltenbreite wie die Zahl, sonst
        # verrutscht die Zeile bei Molekuelen ohne Halogen.
        sig = (f"{'-':>12}" if r.get("sigma_max") is None
               else f"{r['sigma_max']:+12.4f}")
        on = (r.get("vmax_atom") or "?")
        # Bei mehreren Halogenen ist die sigma-Loch-Spalte der GROESSTE Wert.
        # Damit das nicht so aussieht, als gaebe es nur eines, wird hier
        # angehaengt, auf welchem Atom er sitzt und wie viele es insgesamt sind.
        note = ""
        if r.get("n_halogens", 0) > 1:
            note = (f"   {ansi.atom_label(r.get('halogen_atom') or '?')}"
                    f" of {r['n_halogens']}")
        print(f"{r['prefix']:<20}{r['vmin']:>+10.4f}{r['vmax']:>+10.4f}"
              f"{'':>{max(0, 7 - len(on))}}{ansi.atom_label(on)}"
              f"{sig}{r['esp_range']:>9.4f}{note}")
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
