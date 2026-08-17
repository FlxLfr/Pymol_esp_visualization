#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xyzToCube.py
============

Konvertiert Turbomole-``pointval``-Gitterdateien (td.xyz, tp.xyz, ...) in das
Gaussian-Cube-Format, damit sie in PyMOL, VMD, ChimeraX, Avogadro oder Multiwfn
geladen werden koennen.

Hintergrund
-----------
Turbomole schreibt Volumendaten als reine ASCII-Punktwolke::

    #origin           0.000000      0.000000      0.000000
    #vector1          1.000000      0.000000      0.000000
    #vector2          0.000000      1.000000      0.000000
    #vector3          0.000000      0.000000      1.000000
    #grid1  start  -15.000000  delta    0.120000  points    251
    #grid2  start  -15.000000  delta    0.120000  points    251
    #grid3  start  -15.000000  delta    0.120000  points    251
    #title for this grid 111
    #electrostatic potential
    #plotdata
    # cartesian coordinates x,y,z and f(x,y,z)
          -15.00000000   -15.00000000   -15.00000000   -0.00054019
          ...

Jede Zeile enthaelt die vollen Koordinaten -> bei 251^3 Punkten sind das 1.25 GB.
Eine Cube-Datei speichert dieselbe Information mit implizitem Gitter (~200 MB).

Zwei Stolpersteine, die dieses Skript abfaengt:

1. **Achsenreihenfolge.** In der Turbomole-Datei laeuft *x* am schnellsten,
   im Cube-Format laeuft *z* am schnellsten. Ohne Umsortierung erhaelt man ein
   transponiertes, gespiegeltes Molekuel.
2. **Einheiten.** Das Gitter steht in Bohr (atomare Einheiten), die
   Strukturdatei ueblicherweise in Angstrom. Das Skript rechnet die Atome
   standardmaessig um (``--struct-unit angstrom``).

Benutzung
---------
Typischer Aufruf fuer das Brombenzol-Beispiel::

    python xyzToCube.py --struct brombenzol_aro_opti.xyz td.xyz tp.xyz --pymol

Ergebnis: ``td.cube``, ``tp.cube`` und ``esp.pml`` (fertiges PyMOL-Skript).

Wenn PyMOL mit dem vollen 251^3-Gitter zu langsam wird, jeden zweiten Punkt
verwenden::

    python xyzToCube.py --struct brombenzol_aro_opti.xyz td.xyz tp.xyz --stride 2

Als Strukturdatei werden ``.xyz``, ``.mol``, ``.sdf`` und ``.pdb`` akzeptiert -
dieselben Formate wie in render_esp.py. Empfehlung: beiden Skripten *dieselbe*
Datei geben. Sonst stammen die Atome im Cube-Header aus der einen und die
Staebchen in PyMOL aus der anderen Quelle, und eine Abweichung zwischen beiden
faellt nicht auf, weil keine Fehlermeldung kommt.

Nur numpy wird benoetigt (``pip install numpy``).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time

import numpy as np

# ----------------------------------------------------------------------------
# Konstanten
# ----------------------------------------------------------------------------

BOHR_PER_ANGSTROM = 1.8897259886          # CODATA
ANGSTROM_PER_BOHR = 1.0 / BOHR_PER_ANGSTROM

# Ordnungszahlen fuer die gaengigen Elemente. Erweitern falls noetig.
ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
]
SYMBOL_TO_Z = {sym.upper(): i + 1 for i, sym in enumerate(ELEMENTS)}

CHUNK_BYTES = 1 << 24                      # 16 MB Lesepuffer


# ----------------------------------------------------------------------------
# Strukturdatei einlesen
# ----------------------------------------------------------------------------

def _symbol_to_z(sym: str, path: str) -> int:
    """Elementsymbol (oder Ordnungszahl) -> Ordnungszahl."""
    key = sym.strip().capitalize().upper()
    if key in SYMBOL_TO_Z:
        return SYMBOL_TO_Z[key]
    if re.fullmatch(r"\d+", sym.strip()):
        return int(sym)                                # Ordnungszahl statt Symbol
    raise ValueError(
        f"Unbekanntes Element '{sym}' in {path}. "
        f"Bitte die Liste ELEMENTS im Skript ergaenzen."
    )


def _read_xyz(lines, path):
    """xyz-Format: Zeilen ``Symbol x y z``.

    Akzeptiert sowohl das Standardformat (Atomanzahl + Kommentarzeile + Atome)
    als auch eine nackte Koordinatenliste, wie Turbomole-Workflows sie haeufig
    produzieren.
    """
    start = 0
    first = lines[0].strip() if lines else ""
    if re.fullmatch(r"\d+", first):
        start = 2                                     # Anzahl + Kommentar ueberspringen

    atoms = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
        except ValueError:
            continue                                   # Kommentar-/Muellzeile
        atoms.append((_symbol_to_z(parts[0], path), x, y, z))
    return atoms


def _read_molfile(lines, path):
    """MDL-Molfile / SD-File (.mol, .sdf), V2000 und V3000.

    Aufbau V2000::

        Zeile 1   Titel
        Zeile 2   Programmzeile
        Zeile 3   Kommentar
        Zeile 4   Zaehlzeile:  " 12 12  0 ... V2000"
        dann      je Atom:  x  y  z  Symbol  ...      <- Koordinaten ZUERST
        dann      Bindungsblock

    Achtung: die Spaltenreihenfolge ist genau umgekehrt zu xyz. Genau daran
    scheiterte diese Funktion frueher mit "Unbekanntes Element 0.0000".

    Bei SD-Files wird nur der erste Datensatz gelesen (bis ``$$$$``).
    Molfile-Koordinaten sind per Definition in Angstrom.
    """
    atoms = []

    # --- V3000 ---------------------------------------------------------
    if any("V3000" in ln for ln in lines[:8]):
        inside = False
        for line in lines:
            s = line.strip()
            if s.startswith("M  V30 BEGIN ATOM"):
                inside = True
                continue
            if s.startswith("M  V30 END ATOM"):
                break
            if inside and s.startswith("M  V30"):
                # M  V30 <index> <symbol> <x> <y> <z> <aamap> ...
                parts = s.split()
                if len(parts) >= 7:
                    atoms.append((_symbol_to_z(parts[3], path),
                                  float(parts[4]), float(parts[5]),
                                  float(parts[6])))
        return atoms

    # --- V2000 ---------------------------------------------------------
    if len(lines) < 5:
        raise ValueError(f"{path}: zu kurz fuer ein Molfile.")

    counts = lines[3]
    try:
        natoms = int(counts[0:3])
    except ValueError:
        raise ValueError(
            f"{path}: Zaehlzeile (Zeile 4) nicht lesbar: {counts.strip()!r}")

    for line in lines[4:4 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: Atomzeile unvollstaendig: {line.strip()!r}")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        atoms.append((_symbol_to_z(parts[3], path), x, y, z))

    return atoms


def _read_pdb(lines, path):
    """PDB: ATOM-/HETATM-Zeilen, Koordinaten in Angstrom.

    Das Elementsymbol steht in den Spalten 77-78; fehlt es, wird es aus dem
    Atomnamen (Spalten 13-16) abgeleitet.
    """
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        sym = line[76:78].strip()
        if not sym:
            name = line[12:16].strip()
            sym = re.sub(r"[^A-Za-z]", "", name)[:2]
            if len(sym) == 2 and sym.capitalize().upper() not in SYMBOL_TO_Z:
                sym = sym[0]
        atoms.append((_symbol_to_z(sym, path), x, y, z))
    return atoms


def read_structure(path: str, unit: str = "angstrom"):
    """Liest eine Strukturdatei und liefert die Atome in **Bohr**.

    Unterstuetzt:

    ==========  =====================================================
    ``.xyz``    ``Symbol x y z``, mit oder ohne Kopfzeilen
    ``.mol``    MDL-Molfile V2000/V3000 (Koordinaten *vor* dem Symbol)
    ``.sdf``    SD-File, erster Datensatz
    ``.pdb``    ATOM-/HETATM-Zeilen
    ==========  =====================================================

    Rueckgabe: Liste von ``(Z, x, y, z)``.

    ``unit`` gilt nur fuer xyz-Dateien - Molfile und PDB sind per Definition
    in Angstrom, dort wird die Angabe ignoriert.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    if not lines:
        raise ValueError(f"{path} ist leer.")

    ext = os.path.splitext(path)[1].lower()

    if ext in (".mol", ".sdf", ".sd"):
        # SD-File: nur der erste Datensatz
        for i, ln in enumerate(lines):
            if ln.startswith("$$$$"):
                lines = lines[:i]
                break
        atoms = _read_molfile(lines, path)
        file_unit = "angstrom"
    elif ext in (".pdb", ".ent"):
        atoms = _read_pdb(lines, path)
        file_unit = "angstrom"
    else:
        atoms = _read_xyz(lines, path)
        file_unit = unit

    if not atoms:
        raise ValueError(f"Keine Atome in {path} gefunden.")

    if file_unit == "angstrom":
        atoms = [(z, x * BOHR_PER_ANGSTROM,
                     y * BOHR_PER_ANGSTROM,
                     zz * BOHR_PER_ANGSTROM) for (z, x, y, zz) in atoms]
    elif file_unit != "bohr":
        raise ValueError("unit muss 'angstrom' oder 'bohr' sein")

    return atoms


# ----------------------------------------------------------------------------
# Turbomole-Gitterdatei einlesen
# ----------------------------------------------------------------------------

def parse_header(fh):
    """Liest den ``#``-Kopf einer pointval-Datei.

    Rueckgabe: (info-dict, erste_datenzeile). Die erste Datenzeile wurde bereits
    aus dem Stream gelesen und muss vom Aufrufer mitverarbeitet werden.
    """
    info = {
        "origin": np.zeros(3),
        "vectors": np.eye(3),
        "grid": [None, None, None],       # je (start, delta, points)
        "title": "",
        "quantity": "",
    }
    first_data_line = None

    while True:
        line = fh.readline()
        if not line:
            raise ValueError("Datei endet im Header - keine Daten gefunden.")
        if not line.startswith("#"):
            first_data_line = line
            break

        body = line[1:].strip()
        low = body.lower()

        if low.startswith("origin"):
            info["origin"] = np.array([float(v) for v in body.split()[1:4]])
        elif low.startswith("vector"):
            idx = int(body[6]) - 1
            info["vectors"][idx] = np.array([float(v) for v in body.split()[1:4]])
        elif low.startswith("grid"):
            idx = int(body[4]) - 1
            m = re.search(
                r"start\s+(\S+)\s+delta\s+(\S+)\s+points\s+(\d+)", body, re.I)
            if not m:
                raise ValueError(f"Gitterzeile nicht lesbar: {line!r}")
            info["grid"][idx] = (float(m.group(1)),
                                 float(m.group(2)),
                                 int(m.group(3)))
        elif low.startswith("title"):
            info["title"] = body
        elif low in ("density", "electrostatic potential", "plotdata") \
                or "potential" in low or "density" in low:
            if low != "plotdata" and not low.startswith("cartesian"):
                info["quantity"] = body

    if any(g is None for g in info["grid"]):
        raise ValueError("Header unvollstaendig: #grid1/#grid2/#grid3 fehlen.")

    return info, first_data_line


def read_values(path, verbose=True):
    """Liest die 4. Spalte einer pointval-Datei als float32-Array.

    Liest in grossen Bloecken statt Zeile fuer Zeile - fuer die 1.25-GB-Dateien
    ist das etwa eine Groessenordnung schneller.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        info, first_line = parse_header(fh)

        n1 = info["grid"][0][2]
        n2 = info["grid"][1][2]
        n3 = info["grid"][2][2]
        total = n1 * n2 * n3

        values = np.empty(total, dtype=np.float32)
        filled = 0
        t0 = time.time()

        remainder = first_line
        while True:
            chunk = fh.read(CHUNK_BYTES)
            if not chunk:
                break
            chunk = remainder + chunk
            cut = chunk.rfind("\n")
            if cut == -1:                      # extrem lange Zeile - weiterlesen
                remainder = chunk
                continue
            remainder = chunk[cut + 1:]
            block = chunk[:cut]

            tokens = block.split()
            if not tokens:
                continue
            if len(tokens) % 4 != 0:
                raise ValueError(
                    f"{path}: erwartete 4 Spalten pro Zeile, "
                    f"gefunden {len(tokens)} Werte in einem Block."
                )
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            n = arr.size
            if filled + n > total:
                raise ValueError(
                    f"{path}: mehr Datenpunkte als der Header angibt "
                    f"({filled + n} > {total})."
                )
            values[filled:filled + n] = arr
            filled += n

            if verbose:
                pct = 100.0 * filled / total
                sys.stdout.write(f"\r    lese {os.path.basename(path)}: "
                                 f"{pct:5.1f} %")
                sys.stdout.flush()

        # letzte, unvollstaendig gepufferte Zeile
        tokens = remainder.split()
        if tokens:
            if len(tokens) % 4 != 0:
                raise ValueError(f"{path}: letzte Zeile unvollstaendig.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            values[filled:filled + arr.size] = arr
            filled += arr.size

    if verbose:
        sys.stdout.write(f"\r    lese {os.path.basename(path)}: 100.0 %  "
                         f"({filled:,} Punkte in {time.time() - t0:.1f} s)\n")

    if filled != total:
        raise ValueError(
            f"{path}: {filled} Werte gelesen, laut Header erwartet {total}."
        )

    # Turbomole: x laeuft am schnellsten -> Speicherlayout ist [i3, i2, i1]
    data = values.reshape(n3, n2, n1)
    # Cube: z laeuft am schnellsten -> wir wollen [i1, i2, i3]
    data = np.ascontiguousarray(np.transpose(data, (2, 1, 0)))

    return info, data


# ----------------------------------------------------------------------------
# Cube schreiben
# ----------------------------------------------------------------------------

def write_cube(path, info, data, atoms, stride=1, comment=""):
    """Schreibt ein Gaussian-Cube. Alle Laengen in Bohr."""
    if stride > 1:
        data = data[::stride, ::stride, ::stride]

    n = data.shape
    starts = np.array([info["grid"][i][0] for i in range(3)])
    deltas = np.array([info["grid"][i][1] for i in range(3)])
    vecs = info["vectors"]

    # Ursprung des ersten Voxels im kartesischen Raum
    origin = info["origin"] + sum(starts[i] * vecs[i] for i in range(3))
    # Voxelvektoren (mit Stride skaliert)
    voxel = np.array([deltas[i] * stride * vecs[i] for i in range(3)])

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{comment or 'Cube erzeugt mit xyzToCube.py'}\n")
        fh.write(f"{info.get('quantity', '') or 'volumetric data'} | "
                 f"{info.get('title', '')} | Einheiten: Bohr\n")
        fh.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} "
                 f"{origin[2]:12.6f}\n")
        for i in range(3):
            fh.write(f"{n[i]:5d} {voxel[i][0]:12.6f} {voxel[i][1]:12.6f} "
                     f"{voxel[i][2]:12.6f}\n")
        for (znum, x, y, z) in atoms:
            fh.write(f"{znum:5d} {float(znum):12.6f} {x:12.6f} {y:12.6f} "
                     f"{z:12.6f}\n")

        # Werte: z am schnellsten, 6 pro Zeile.
        # Ein vorkompiliertes Formatmuster pro z-Reihe ist deutlich schneller
        # als eine Schleife ueber alle 15.8 Mio. Einzelwerte.
        nz = n[2]
        n_full, rest = divmod(nz, 6)
        row_fmt = ("%13.5E" * 6 + "\n") * n_full
        if rest:
            row_fmt += "%13.5E" * rest + "\n"

        flat = np.ascontiguousarray(data).reshape(n[0] * n[1], nz)
        buf = []
        for idx in range(flat.shape[0]):
            buf.append(row_fmt % tuple(flat[idx].tolist()))
            if len(buf) >= 4096:                    # gebuendelt schreiben
                fh.write("".join(buf))
                buf.clear()
        if buf:
            fh.write("".join(buf))

    return n, origin, voxel


# ----------------------------------------------------------------------------
# PyMOL-Skript erzeugen
# ----------------------------------------------------------------------------

PML_TEMPLATE = """# --------------------------------------------------------------
# esp.pml - ESP auf Elektronendichte-Isoflaeche
# Erzeugt von xyzToCube.py
# Start:  pymol esp.pml       (oder in PyMOL:  @esp.pml)
# --------------------------------------------------------------

reinitialize

# 1) Struktur und Volumendaten laden
load {struct}, mol
{load_density}
load {esp_cube}, esp

# 2) Molekuel als Staebchen
hide everything
show sticks, mol
set stick_radius, 0.12
color grey70, mol and elem C
util.cnc mol

# 3) Isoflaeche der Elektronendichte bei rho = {iso} a.u.
#    (Politzer/Murray-Konvention fuer die "Molekueloberflaeche")
{isosurface}

# 4) Farbrampe fuer das ESP; Werte in Hartree/e (a.u.)
#    {vmin} .. {vmax} a.u.  entspricht {kvmin:.0f} .. {kvmax:.0f} kJ/(mol*e)
ramp_new espramp, esp, [{vmin}, 0.0, {vmax}], [red, white, blue]

# 5) ESP auf die Oberflaeche mappen
set surface_color, espramp, {surface_target}
set surface_quality, 1

#    Transparenz: 0 = opak (kraftigste Farben, Staebchen unsichtbar),
#    0.15 = Standard (Molekuelgeruest scheint durch),
#    ab ca. 0.3 wird es unleserlich, weil man durch das ganze Molekuel schaut.
set transparency, {transparency}
set transparency_mode, 2
set two_sided_lighting, on

# 6) Darstellung / Rendering
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 0
set specular, 0.2
set ambient, 0.15
orient mol
zoom mol, 2.0

# 7) Hochaufloesendes Bild
# ray 2400, 1800
# png esp.png, dpi=300
"""


def write_pymol_script(path, struct, density_cube, esp_cube, vmin, vmax,
                       iso=0.001, transparency=0.15):
    if density_cube:
        load_density = f"load {density_cube}, dens"
        isosurface = f"isosurface surf, dens, {iso}"
        surface_target = "surf"
    else:
        load_density = "# keine Dichtedatei vorhanden"
        isosurface = ("# Ersatz: van-der-Waals-Oberflaeche aus der Struktur.\n"
                      "# Fuer die Politzer/Murray-Konvention (rho = 0.001 a.u.)\n"
                      "# waere die Dichtedatei td.cube noetig.\n"
                      "show surface, mol\n"
                      "set surface_solvent, 0")
        surface_target = "mol"

    # 1 Hartree/e = 2625.5 kJ/(mol*e)
    text = PML_TEMPLATE.format(
        struct=struct,
        load_density=load_density,
        esp_cube=esp_cube,
        isosurface=isosurface,
        surface_target=surface_target,
        iso=iso,
        transparency=transparency,
        vmin=vmin, vmax=vmax,
        kvmin=vmin * 2625.4996, kvmax=vmax * 2625.4996,
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


# ----------------------------------------------------------------------------
# Hauptprogramm
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Turbomole-pointval-Gitterdateien (td.xyz, tp.xyz) "
                    "nach Gaussian-Cube konvertieren.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Beispiel:\n"
               "  python xyzToCube.py --struct brombenzol_aro_opti.xyz "
               "td.xyz tp.xyz --pymol\n",
    )
    p.add_argument("grids", nargs="+",
                   help="Turbomole-Gitterdateien (z.B. td.xyz tp.xyz)")
    p.add_argument("--struct", "-s", required=True,
                   help="Strukturdatei: .xyz, .mol, .sdf oder .pdb")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom",
                   help="Einheit der Strukturdatei (Standard: angstrom); gilt nur fuer .xyz - .mol/.sdf/.pdb sind immer Angstrom")
    p.add_argument("--outdir", "-o", default=None,
                   help="Ausgabeverzeichnis (Standard: neben der Eingabe)")
    p.add_argument("--stride", type=int, default=1,
                   help="Nur jeden n-ten Gitterpunkt schreiben "
                        "(2 => 8x kleinere Datei, Standard: 1)")
    p.add_argument("--pymol", action="store_true",
                   help="Zusaetzlich ein fertiges esp.pml schreiben")
    p.add_argument("--esp-range", type=float, default=0.03,
                   help="Halbe Breite der ESP-Farbskala in a.u. "
                        "(Standard: 0.03 => -0.03 .. +0.03)")
    p.add_argument("--iso", type=float, default=0.001,
                   help="Isowert der Dichteflaeche in a.u. (Standard: 0.001)")
    p.add_argument("--transparency", type=float, default=0.15,
                   help="Oberflaechentransparenz im esp.pml, 0..1 "
                        "(Standard: 0.15; 0 = opak)")
    p.add_argument("--quiet", "-q", action="store_true")
    args = p.parse_args(argv)

    verbose = not args.quiet

    if verbose:
        print("=" * 70)
        print("xyzToCube.py - Turbomole pointval  ->  Gaussian Cube")
        print("=" * 70)

    atoms = read_structure(args.struct, unit=args.struct_unit)
    if verbose:
        print(f"[1] Struktur: {args.struct} -> {len(atoms)} Atome "
              f"(eingelesen als {args.struct_unit}, gespeichert als Bohr)")

    written = {}
    for gpath in args.grids:
        if verbose:
            print(f"[2] Gitterdatei: {gpath}")
        info, data = read_values(gpath, verbose=verbose)

        n1, n2, n3 = (info["grid"][i][2] for i in range(3))
        if verbose:
            print(f"    Gitter {n1} x {n2} x {n3}, "
                  f"delta = {info['grid'][0][1]} Bohr, "
                  f"Groesse = '{info['quantity'] or 'unbekannt'}'")
            print(f"    Wertebereich: {data.min():+.6g} .. {data.max():+.6g}")

        base = os.path.splitext(os.path.basename(gpath))[0]
        outdir = args.outdir or os.path.dirname(os.path.abspath(gpath))
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, base + ".cube")

        shape, origin, voxel = write_cube(
            outpath, info, data, atoms, stride=args.stride,
            comment=f"{info['quantity'] or base} - konvertiert aus {os.path.basename(gpath)}",
        )
        if verbose:
            mb = os.path.getsize(outpath) / 1024 ** 2
            print(f"    -> {outpath}  ({shape[0]}x{shape[1]}x{shape[2]}, "
                  f"{mb:.1f} MB)")

        q = (info["quantity"] or "").lower()
        if "potential" in q:
            written["esp"] = outpath
        elif "density" in q:
            written["density"] = outpath
        else:
            written.setdefault("other", []).append(outpath)

    if args.pymol:
        outdir = args.outdir or os.path.dirname(os.path.abspath(args.grids[0]))
        pml = os.path.join(outdir, "esp.pml")
        esp_cube = written.get("esp")
        if esp_cube is None:
            print("    ! Keine ESP-Datei erkannt - esp.pml wird uebersprungen.",
                  file=sys.stderr)
        else:
            write_pymol_script(
                pml,
                struct=os.path.relpath(os.path.abspath(args.struct), outdir),
                density_cube=(os.path.basename(written["density"])
                              if "density" in written else None),
                esp_cube=os.path.basename(esp_cube),
                vmin=-args.esp_range, vmax=args.esp_range, iso=args.iso,
                transparency=args.transparency,
            )
            if verbose:
                print(f"[3] PyMOL-Skript: {pml}")
                print(f"    Start mit:  pymol {os.path.basename(pml)}")

    if verbose:
        print("Fertig.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
