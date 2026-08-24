#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_esp.py
=============

Erzeugt einen standardisierten Satz ESP-Bilder aus einem Paar Cube-Dateien
(Elektronendichte + elektrostatisches Potential), vollautomatisch und ohne
einen einzigen Mausklick.

Das ESP wird auf die Isoflaeche der Elektronendichte bei rho = 0.001 a.u.
abgebildet (Konvention nach Politzer/Murray) und aus drei fest definierten
Blickrichtungen gerendert:

    pi      senkrecht auf die Molekuelebene   -> zeigt das pi-System
    edge    in der Molekuelebene              -> Profil, C-X-Achse waagerecht
    sigma   entlang der C-X-Achse von aussen  -> zeigt das sigma-Loch frontal

Die Orientierung wird aus der Geometrie berechnet (Traegheitsachsen +
Kohlenstoff-Halogen-Achse), NICHT ueber PyMOLs ``orient``. Dadurch liefern
verschiedene Molekuele reproduzierbar dieselbe Ausrichtung.


Aufruf
------
    pymol -cq render_esp.py -- --density td.cube --esp tp.cube \\
                               --struct brombenzol_aro_opti.mol \\
                               --prefix brombenzol

Ohne Argumente sucht das Skript im aktuellen Ordner nach ``td.cube``,
``tp.cube`` und einer Struktur (.mol/.sdf/.xyz).


Farbskala
---------
Standardmaessig wird der ESP-Bereich aus den Daten *auf der Isoflaeche*
bestimmt und symmetrisch auf einen glatten Wert aufgerundet. Der verwendete
Wert wird ins Log und in die Datei ``<prefix>_settings.txt`` geschrieben.

!! Fuer den direkten Vergleich mehrerer Molekuele muss die Skala fest sein.
   Dazu einmal alle Molekuele mit --esp-range auto durchlaufen lassen, den
   groessten gemeldeten Wert notieren und danach alle erneut mit
   --esp-range <wert> rendern.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

import ansi
import xyzToCube                    # Elementliste (siehe Z_SYMBOL) und Schale
# Die Schalenauswertung (rho = iso) und die Farbskala stehen in xyzToCube.py,
# weil das Konvertierskript sie fuer sein eigenes esp.pml ebenfalls braucht.
# Beide Skripte benutzen damit dieselbe Definition von "auf der Isoflaeche".
from xyzToCube import esp_statistics, nice_range, shell_mask, SHELL_TOL_FACTOR
from constants import BOHR_PER_ANGSTROM, HARTREE_TO_KCAL, HARTREE_TO_KJ


# ----------------------------------------------------------------------------
# Cube einlesen (nur fuer die Statistik; PyMOL laedt die Dateien selbst)
# ----------------------------------------------------------------------------

def read_cube(path):
    """Liest ein Gaussian-Cube. Rueckgabe: (werte3d, atome, origin, voxel)."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        fh.readline()
        parts = fh.readline().split()
        natoms = int(parts[0])
        origin = np.array([float(v) for v in parts[1:4]])

        n = []
        voxel = []
        for _ in range(3):
            p = fh.readline().split()
            n.append(int(p[0]))
            voxel.append([float(v) for v in p[1:4]])
        n = np.array(n)
        voxel = np.array(voxel)

        atoms = []
        for _ in range(natoms):
            p = fh.readline().split()
            atoms.append((int(p[0]),
                          float(p[2]), float(p[3]), float(p[4])))

        values = np.fromstring(fh.read(), sep=" ", dtype=np.float64)

    if values.size != int(np.prod(n)):
        raise ValueError(f"{path}: erwartet {int(np.prod(n))} Werte, "
                         f"gelesen {values.size}")
    return values.reshape(n), atoms, origin, voxel


def shell_points(density, esp, origin, voxel, iso=0.001,
                 tol_factor=SHELL_TOL_FACTOR):
    """Koordinaten und ESP-Werte der Gitterpunkte auf der rho=iso-Schale.
    Nur die Schalenpunkte werden materialisiert, nicht das ganze Gitter -
    bei 251^3 waere ein volles Koordinatenfeld sonst mehrere hundert MB.
    """
    mask = shell_mask(density, iso, tol_factor)
    idx = np.argwhere(mask)                       # (N, 3) Gitterindizes
    pos = origin + idx @ voxel                    # (N, 3) kartesisch, Bohr
    return pos, esp[mask]


def halogen_axes(atoms):
    """Alle Halogene mit ihrer C-X-Achse.

    Ein Molekuel kann mehr als ein Halogen tragen - Jedes wird einzeln ausgewertet.

    Rueckgabe: Liste von dicts mit
      index   0-basierter Atomindex des Halogens
      symbol  Elementsymbol
      label   Symbol + 1-basierte Nummer, z.B. "Cl12" (wie in der Ausgabe)
      pos     Koordinaten des Halogens
      axis    normierte C->X-Achse (zeigt zum sigma-Loch)
      r_limit maximaler Abstand, in dem die eigene Oberflaeche des Halogens
              liegen kann (Bohr) - siehe unten
    Halogene ohne gebundenen Kohlenstoff in Reichweite werden uebersprungen.

    Warum ``r_limit``: bei gefalteten Molekuelen zeigt der Kegel um die C-X-
    Achse nicht ins Leere, sondern auf einen anderen Molekuelteil. Ohne
    Abstandsgrenze wird dessen Oberflaeche mitgemessen. Triazolam ist genau so
    ein Fall - der Kegel um Cl21 trifft die Methylgruppe am Triazolring, und
    die lieferte ein "sigma-Loch" von +18.8 statt +10.5 kcal/(mol*e). Die
    rho=0.001-Flaeche eines Halogens liegt bei etwa 1.1 bis 1.2 vdW-Radien;
    der Faktor 1.6 laesst reichlich Luft und schliesst alles Weitere aus.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])
    carbons = [i for i, z in enumerate(znums) if z == 6]
    out = []
    if not carbons:
        return out
    for hi, z in enumerate(znums):
        if int(z) not in HALOGENS:
            continue
        d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
        ci = carbons[int(np.argmin(d))]
        axis = coords[hi] - coords[ci]
        n = np.linalg.norm(axis)
        if n < 1e-6:
            continue
        out.append({"index": hi,
                    "symbol": HALOGENS[int(z)],
                    "label": f"{HALOGENS[int(z)]}{hi + 1}",
                    "pos": coords[hi],
                    "axis": axis / n,
                    "r_limit": (1.6 * VDW_ANGSTROM[int(z)]
                                * BOHR_PER_ANGSTROM)})
    return out


def local_extrema(pos, vals, atoms,
                  cone_cos=0.80, belt_cos=0.35, belt_factor=1.5):
    """Regionsaufgeloeste Oberflaechen-Extrema.

    Warum das noetig ist: das *globale* V_S,max einer Aryl-Halogenid-Oberflaeche
    sitzt fast immer auf den Ring-Wasserstoffen, nicht auf dem Halogen. Bei
    Brombenzol und Iodbenzol liefert der globale Wert deshalb zweimal dasselbe
    C-H (+0.031 a.u.), waehrend sich die sigma-Loecher, um die es eigentlich
    geht, um fast den Faktor zwei unterscheiden. Fuer jeden Vergleich von
    Halogenbruecken-Donoren braucht man das *lokale* Maximum auf dem Halogen.

    Regionen:
      sigma   Kappe um die verlaengerte C-X-Achse (Oeffnungswinkel aus cone_cos)
      belt    Guertel senkrecht dazu, in Halogennaehe

    Rueckgabe: dict mit den Kennwerten in atomaren Einheiten; die
    halogenbezogenen Eintraege fehlen, wenn das Molekuel kein Halogen enthaelt.
    """
    out = {}
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])

    # --- Region des globalen Extremums benennen ------------------------
    for tag, i in (("vmax", int(np.argmax(vals))), ("vmin", int(np.argmin(vals)))):
        d = np.linalg.norm(coords - pos[i], axis=1)
        j = int(np.argmin(d))
        out[f"{tag}_atom"] = f"{z_symbol(int(znums[j]))}{j + 1}"

    # --- Halogenregionen, eines nach dem anderen -------------------------
    entries = []
    for hal in halogen_axes(atoms):
        e = {k: hal[k] for k in ("index", "symbol", "label")}
        axis = hal["axis"]

        rel = pos - hal["pos"]
        r = np.linalg.norm(rel, axis=1)
        r[r == 0] = 1e-9
        cos = (rel @ axis) / r

        # Abstandsgrenze: sonst zaehlen bei gefalteten Molekuelen Punkte auf
        # einem ganz anderen Molekuelteil zur Kappe (siehe halogen_axes).
        cap = (cos > cone_cos) & (r < hal["r_limit"])
        if cap.sum() >= 5:
            e["sigma_max"] = float(vals[cap].max())
            e["sigma_points"] = int(cap.sum())
            r_cap = float(r[cap].mean())
            belt = (np.abs(cos) < belt_cos) & (r < belt_factor * r_cap)
            if belt.sum() >= 5:
                e["belt_min"] = float(vals[belt].min())
                e["belt_points"] = int(belt.sum())
        entries.append(e)

    if entries:
        out["halogens"] = entries
    return out


def rank_halogens(entries):
    """Halogeneintraege nach sigma-Loch absteigend sortieren.

    Eintraege ohne auswertbares sigma-Loch wandern ans Ende. Die Reihenfolge
    bestimmt die Ausgabe und - ueber den ersten Eintrag - die Orientierung der
    sigma-Ansicht.
    """
    return sorted(entries,
                  key=lambda e: -e.get("sigma_max", -np.inf))


def promote_primary(out):
    """Kennwerte des staerksten sigma-Lochs zusaetzlich flach in ``out`` legen.

    Damit bleiben Aufrufer, die nur EINEN Wert erwarten (CSV-Spalte
    ``sigma_hole_au``, Zusammenfassungstabelle), unveraendert lauffaehig; bei
    einem einzelnen Halogen ist das Ergebnis bitgleich zu vorher.
    """
    entries = out.get("halogens") or []
    if not entries:
        return out
    entries = rank_halogens(entries)
    out["halogens"] = entries
    first = entries[0]
    out["halogen"] = first["symbol"]
    out["halogen_atom"] = first["label"]
    out["halogen_index"] = first["index"]
    for k in ("sigma_max", "sigma_points",
              "sigma_angle", "sigma_method", "belt_min", "belt_points"):
        if k in first:
            out[k] = first[k]
    return out



def _trilinear(vol, origin, delta, pts):
    """Trilineare Interpolation auf einem achsparallelen, regelmaessigen Gitter."""
    f = (pts - origin) / delta
    i0 = np.floor(f).astype(int)
    frac = f - i0
    n = np.array(vol.shape)
    i0 = np.clip(i0, 0, n - 2)
    frac = np.clip(f - i0, 0.0, 1.0)

    out = np.zeros(len(pts))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((frac[:, 0] if dx else 1 - frac[:, 0]) *
                     (frac[:, 1] if dy else 1 - frac[:, 1]) *
                     (frac[:, 2] if dz else 1 - frac[:, 2]))
                out += w * vol[i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz]
    return out


def _cone_directions(axis, cone_cos, n=400):
    """Richtungen in einer Kappe um ``axis``; die erste ist die Achse selbst.

    Der Rest ist eine Fibonacci-Spirale auf der Kugelkappe - gleichmaessige
    Belegung ohne Haeufung an der Achse, wie sie bei Kugelkoordinaten auftraete.

    Warum die Achse gesondert: das ``+ 0.5`` in ``k`` ist die Mittelpunktsregel
    fuer FLAECHENGLEICHE Verteilung - jeder Strahl sitzt in der Mitte eines
    gleich grossen Rings. Das ist richtig, wenn man ueber die Kappe mittelt.
    Wir suchen aber ein MAXIMUM, und das sitzt beim sigma-Loch genau auf der
    Achse. Mit dem Versatz war der innerste Strahl 1.281 Grad davon entfernt,
    die Achse selbst wurde nie ausgewertet, und jedes achsensymmetrische
    Molekuel meldete stur "1.3 Grad" - eine Untergrenze des Abtastrasters,
    keine Messung. Der Fehler im Wert war klein (4-Bromacetophenon: 0.008
    kcal/(mol*e)), die Meldung aber irrefuehrend.
    """
    axis = np.asarray(axis, dtype=float)

    # n-1 Spiralrichtungen; die Achse kommt als erste dazu
    m = max(1, n - 1)
    k = np.arange(m) + 0.5
    cosv = 1.0 - (1.0 - cone_cos) * k / m          # cone_cos .. 1
    phi = np.pi * (1 + 5 ** 0.5) * k
    sinv = np.sqrt(np.maximum(0.0, 1 - cosv ** 2))

    # orthonormale Basis um axis
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, axis)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    spiral = (cosv[:, None] * axis
              + (sinv * np.cos(phi))[:, None] * e1
              + (sinv * np.sin(phi))[:, None] * e2)
    return np.vstack([axis[None, :], spiral])


def sigma_hole_interpolated(density, esp, origin, voxel, atoms, iso=0.001,
                            cone_cos=0.80, n_rays=400, dr=0.02, r_max=14.0):
    """sigma-Loch ohne Abhaengigkeit von der Gitterausrichtung.

    Das Problem der punktbasierten Auswertung: das sigma-Loch ist ein Gipfel
    *auf* der C-X-Achse. Ob ein Gitterpunkt dort UND gleichzeitig innerhalb der
    duennen rho=iso-Schale liegt, ist Zufall. Bei Brombenzol lag der beste
    Punkt auf dem 126^3-Gitter 1.14 Bohr neben der Achse - der Wert kam 28 %
    zu niedrig heraus, obwohl 144 Punkte in der Kappe lagen.

    Stattdessen hier: Strahlen vom Halogen aus in eine Kappe um die Achse, auf
    jedem Strahl der Radius gesucht, bei dem rho die Isoflaeche schneidet, und
    dort V ausgewertet - beides trilinear interpoliert. Das Ergebnis haengt
    nicht mehr davon ab, wo die Gitterpunkte zufaellig liegen.

    Traegt das Molekuel mehrere Halogene, wird jedes einzeln abgetastet.

    Rueckgabe: dict {Atomindex: {sigma_max, sigma_angle, sigma_method}};
    leeres dict, wenn kein Halogen auswertbar ist.
    """
    diag = np.diag(voxel)
    if not np.allclose(voxel, np.diag(diag)):
        return {}                         # nicht achsparallel - Fallback
    delta = diag

    result = {}

    for hal in halogen_axes(atoms):
        axis = hal["axis"]
        origin_atom = hal["pos"]
        dirs = _cone_directions(axis, cone_cos, n_rays)
        radii = np.arange(1.0, min(r_max, hal["r_limit"]), dr)

        best_v, best_cos = None, None
        for d in dirs:
            pts = origin_atom + radii[:, None] * d[None, :]
            rho = _trilinear(density, origin, delta, pts)
            # INNERSTER Schnittpunkt, von innen nach aussen: der Strahl startet
            # tief in der Dichte des Halogens, das erste Unterschreiten von iso
            # ist dessen eigene Oberflaeche. Frueher wurde der aeusserste
            # Schnittpunkt genommen; bei gefalteten Molekuelen taucht der
            # Strahl dahinter in einen anderen Molekuelteil ein und misst
            # dessen Oberflaeche.
            if rho[0] < iso:
                continue                  # Strahl startet schon ausserhalb
            below = np.nonzero(rho < iso)[0]
            if below.size == 0:
                continue                  # Flaeche innerhalb r_limit nicht getroffen
            j = below[0] - 1
            if j < 0 or j + 1 >= len(radii):
                continue
            # lineare Interpolation des Radius am Isowert
            r0, r1 = radii[j], radii[j + 1]
            y0, y1 = rho[j], rho[j + 1]
            rs = r0 + (iso - y0) * (r1 - r0) / (y1 - y0) if y1 != y0 else r0
            v = float(_trilinear(esp, origin, delta,
                                 (origin_atom + rs * d)[None, :])[0])
            if best_v is None or v > best_v:
                best_v, best_cos = v, float(np.dot(d, axis))

        if best_v is None:
            continue
        result[hal["index"]] = {
            "sigma_max": best_v,
            "sigma_angle": float(np.degrees(np.arccos(min(1.0, best_cos)))),
            "sigma_method": "interpoliert"}

    return result


# ----------------------------------------------------------------------------
# Orientierung aus der Geometrie
# ----------------------------------------------------------------------------

HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}

# ----------------------------------------------------------------------------
# Farbrampen
#
# Beide behalten dieselbe Konvention: ROT = negativ, BLAU = positiv, und die
# Mitte der Skala ist V = 0. Die Regenbogenrampe schiebt nur Gelb/Gruen/Cyan
# dazwischen, statt ueber Weiss zu gehen. Dadurch bleiben die Bilder beider
# Rampen an den Enden vergleichbar - waere die Regenbogenskala umgedreht
# (blau = negativ, wie in manchen Programmen), koennte man die zwei Saetze
# nicht nebeneinanderlegen.
#
# Nutzen der Regenbogenrampe: rot-weiss-blau hat in der Mitte kaum
# Farbaufloesung, schwach polare Bereiche sehen alle gleich weiss aus. Der
# Regenbogen loest genau dort auf. Preis: er ist nicht perzeptuell gleichmaessig
# und erzeugt Kanten, wo keine sind - fuer eine quantitative Aussage bleibt
# rot-weiss-blau die ehrlichere Darstellung.
# ----------------------------------------------------------------------------

RAMP_PYMOL = {
    "redblue": ["red", "white", "blue"],
    "rainbow": ["red", "yellow", "green", "cyan", "blue"],
}
RAMP_HEX = {
    "redblue": ["#d40000", "#ffffff", "#0030d4"],
    "rainbow": ["#d40000", "#f0e000", "#00a000", "#00c8d4", "#0030d4"],
}


def ramp_levels(rng, rainbow=False):
    """Stuetzstellen und Farben fuer ``cmd.ramp_new``.

    Rueckgabe: (levels, colors) - gleich lang, symmetrisch um 0.
    """
    name = "rainbow" if rainbow else "redblue"
    colors = RAMP_PYMOL[name]
    n = len(colors)
    levels = [-rng + 2.0 * rng * i / (n - 1) for i in range(n)]
    return levels, colors

# van-der-Waals-Radien nach Bondi (J. Phys. Chem. 1964, 68, 441), Angstrom.
# Nur als Groessenordnung fuer die Abstandsgrenze der sigma-Loch-Suche.
VDW_ANGSTROM = {9: 1.47, 17: 1.75, 35: 1.85, 53: 1.98}

# Rueckuebersetzung Ordnungszahl -> Symbol, abgeleitet aus derselben Liste, mit
# der xyzToCube.py in die Gegenrichtung uebersetzt. Frueher stand hier eine
# zweite, handgepflegte Tabelle mit 20 Eintraegen: ein Molekuel mit einem
# Element, das dort fehlte, wurde sauber konvertiert und anschliessend als
# "Z13" statt "Al13" beschriftet. Eine Quelle, beide Richtungen.
Z_SYMBOL = {i + 1: sym for i, sym in enumerate(xyzToCube.ELEMENTS)}


def z_symbol(z):
    return Z_SYMBOL.get(int(z), f"Z{int(z)}")


def molecular_frame(atoms, halogen_index=None):
    """Bestimmt ein reproduzierbares Molekuelkoordinatensystem.

    ``halogen_index`` waehlt bei mehreren Halogenen aus, welches die Achse
    festlegt - und damit, auf welches sigma-Loch die sigma-Ansicht blickt.
    Ohne Angabe wird das erste Halogen in der Atomliste genommen; render_all
    uebergibt das Halogen mit dem staerksten sigma-Loch.

    Rueckgabe: (normal, axis, sigma_axis, center)
      normal      Flaechennormale (kleinste Traegheitsausdehnung, Schweratome)
      axis        IN DIE EBENE PROJIZIERTE C->Halogen-Achse. Sie spannt mit
                  ``normal`` ein sauberes Rechtssystem auf und richtet die
                  pi- und edge-Ansicht aus.
      sigma_axis  die ECHTE C->Halogen-Achse, unprojiziert
      center      geometrischer Mittelpunkt aller Atome
    Alle Vektoren normiert, Koordinaten in denselben Einheiten wie ``atoms``.

    Warum zwei Achsen: bei planaren Molekuelen sind beide identisch, die
    Projektion ist dort reine Rundungskosmetik. Bei nicht-planaren Molekuelen
    dreht sie die Achse aber tatsaechlich weg - bei Triazolam steht die
    C-Cl21-Bindung 42.9 Grad aus der Ausgleichsebene heraus, und die
    sigma-Ansicht blickte entsprechend 42.9 Grad am sigma-Loch vorbei. Fuer die
    sigma-Ansicht ist deshalb ``sigma_axis`` zu verwenden.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])
    center = coords.mean(axis=0)

    heavy = coords[znums > 1]
    if len(heavy) < 3:
        heavy = coords
    centered = heavy - heavy.mean(axis=0)

    # Hauptachsen ueber Singulaerwertzerlegung
    _, sing, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[2]                                   # kleinste Ausdehnung
    long_axis = vt[0]                                # groesste Ausdehnung

    # C-Halogen-Achse suchen
    axis = None
    hal_idx = [i for i, z in enumerate(znums) if z in HALOGENS]
    if hal_idx:
        hi = halogen_index if halogen_index in hal_idx else hal_idx[0]
        carbons = [i for i, z in enumerate(znums) if z == 6]
        if carbons:
            d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
            ci = carbons[int(np.argmin(d))]
            axis = coords[hi] - coords[ci]           # C -> X, zeigt zum sigma-Loch

    if axis is None:
        axis = long_axis.copy()

    axis = axis / np.linalg.norm(axis)
    normal = normal / np.linalg.norm(normal)
    sigma_axis = axis.copy()                 # unveraendert, fuer die sigma-Ansicht

    # axis in die Ebene legen - nur fuer pi und edge, die ein Rechtssystem
    # mit der Normalen brauchen.
    axis = axis - normal * float(np.dot(axis, normal))
    if np.linalg.norm(axis) < 1e-6:
        axis = long_axis
    axis = axis / np.linalg.norm(axis)

    return normal, axis, sigma_axis, center


def view_matrix(forward, up):
    """Rotationsmatrix fuer PyMOLs set_view.

    ``forward`` zeigt vom Molekuel zur Kamera, ``up`` nach oben im Bild.
    Zeilen der Matrix sind die Kamera-Basisvektoren im Weltsystem.
    """
    z = np.asarray(forward, dtype=float)
    z = z / np.linalg.norm(z)
    up = np.asarray(up, dtype=float)
    up = up - z * float(np.dot(up, z))
    if np.linalg.norm(up) < 1e-8:                    # up parallel zu z
        alt = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(alt, z))) > 0.9:
            alt = np.array([1.0, 0.0, 0.0])
        up = alt - z * float(np.dot(alt, z))
    up = up / np.linalg.norm(up)
    right = np.cross(up, z)
    right = right / np.linalg.norm(right)
    # PyMOL erwartet in set_view die Matrix, deren SPALTEN die
    # Kamera-Basisvektoren im Weltsystem sind - also die Transponierte
    # der Zeilenform. Empirisch geprueft (siehe SOP, Abschnitt Ansichten).
    return np.array([right, up, z]).T


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------

def ensure_pymol():
    """Liefert ``cmd``; startet PyMOL headless, falls noch nicht gestartet.

    So laeuft das Skript sowohl mit ``python render_esp.py ...`` als auch mit
    ``pymol -cq render_esp.py -- ...``.
    """
    import pymol
    if not getattr(pymol, "_cmd", None) or not hasattr(pymol, "cmd"):
        pymol.finish_launching(["pymol", "-qc"])
    else:
        try:
            pymol.cmd.get_version()
        except Exception:
            pymol.finish_launching(["pymol", "-qc"])
    return pymol.cmd


def render_all(args):
    cmd = ensure_pymol()

    # --- Daten fuer die Statistik ---------------------------------------
    dens, atoms, origin, voxel = read_cube(args.density)
    esp, _, _, _ = read_cube(args.esp)
    if dens.shape != esp.shape:
        raise SystemExit("Dichte- und ESP-Cube haben unterschiedliche Gitter.")

    vmin, vmax, npts = esp_statistics(dens, esp, iso=args.iso)
    if npts == 0:
        raise SystemExit(f"Keine Gitterpunkte bei rho = {args.iso} gefunden. "
                         f"Isowert pruefen.")

    if args.esp_range == "auto":
        rng = nice_range(vmin, vmax)
        how = "automatisch aus den Daten"
    else:
        rng = float(args.esp_range)
        how = "fest vorgegeben"

    # Regionsaufgeloeste Kennwerte: das globale Maximum sitzt bei Arylhalogeniden
    # auf den Ring-Wasserstoffen, nicht auf dem Halogen.
    pos, vals = shell_points(dens, esp, origin, voxel, iso=args.iso)
    loc = local_extrema(pos, vals, atoms)

    # Das sigma-Loch wird bevorzugt ueber Strahlen mit Interpolation bestimmt.
    # Die punktbasierte Variante haengt davon ab, ob zufaellig ein Gitterpunkt
    # nahe der C-X-Achse UND in der duennen Isoschale liegt; bei Brombenzol
    # ergibt sie auf demselben Gitter +7.9 statt +10.1 kcal/(mol*e).
    ray = sigma_hole_interpolated(dens, esp, origin, voxel, atoms, iso=args.iso)
    for e in loc.get("halogens", []):
        if e["index"] in ray:
            e.update(ray[e["index"]])
    promote_primary(loc)
    hals = loc.get("halogens", [])
    spacing = float(np.max(np.abs(np.diag(voxel))))

    def _fmt(v):
        return (f"{v:+.4f} a.u.  = {v*HARTREE_TO_KJ:+7.1f} kJ/(mol*e)"
                f"  = {v*HARTREE_TO_KCAL:+6.1f} kcal/(mol*e)")

    print(f"  ESP auf der rho={args.iso}-Schale ({npts} Punkte):")
    print(f"    V_S,min = {_fmt(vmin)}   auf "
          f"{ansi.atom_label(loc.get('vmin_atom', '?'))}")
    print(f"    V_S,max = {_fmt(vmax)}   auf "
          f"{ansi.atom_label(loc.get('vmax_atom', '?'))}")
    # Pro Halogen ein Block. Bei genau einem Halogen ist die Ausgabe dieselbe
    # wie vorher; ab zwei bekommt jedes seine eigene Zeile, und das Halogen,
    # auf das die sigma-Ansicht blickt, ist markiert.
    for e in hals:
        if "sigma_max" not in e:
            print(f"  Lokal am Halogen ({ansi.element(e['symbol'])}"
                  f" {e['index'] + 1}): kein auswertbares sigma-Loch "
                  f"(zu wenige Oberflaechenpunkte in der Kappe)")
            continue
        head = (f"  Lokal am Halogen ({ansi.element(e['symbol'])}):"
                if len(hals) == 1
                else f"  Lokal an {ansi.atom_label(e['label'])}:")
        if len(hals) > 1 and e is hals[0]:
            head += "   <- Achse der sigma-Ansicht"
        print(head)
        tag = e.get("sigma_method", "punktbasiert")
        extra = (f"   [{tag}, {e['sigma_angle']:.1f} Grad zur Achse]"
                 if "sigma_angle" in e
                 else f"   [{tag}, {e.get('sigma_points', 0)} Punkte]")
        print(f"    sigma-Loch  = {_fmt(e['sigma_max'])}{extra}")
        if "belt_min" in e:
            print(f"    Guertel     = {_fmt(e['belt_min'])}"
                  f"   [{e['belt_points']} Punkte]")
    if hals and any("sigma_max" in e for e in hals):
        # Auch das Strahlverfahren kann die Isoflaeche auf einem groben Gitter
        # nur so genau lokalisieren, wie die Dichte dort aufgeloest ist.
        if spacing > 0.30:
            print(f"    ! Gitterabstand {spacing:.2f} Bohr - fuer einen "
                  f"belastbaren sigma-Loch-Wert zu grob;")
            print(f"      erwartungsgemaess einige Prozent zu niedrig. "
                  f"Feiner rechnen (kleineres --stride).")
        # Hinweis, dass V_S,max nicht auf dem Halogen liegt: bewusst nicht
        # ausgegeben. Bei Arylhalogeniden trifft das praktisch immer zu, die
        # Meldung waere also bei jedem Molekuel identisch und damit wertlos.
        # Die Information steckt bereits in der Ortsangabe hinter V_S,max
        # ("auf H5") und im separat ausgewiesenen sigma-Loch. Erklaerung dazu
        # in docs/ESP_Visualization_Background.docx, Abschnitt 2.1
        # "Which number describes the sigma-hole - Not V_S,max".
    print(f"  Farbskala: +/- {rng:.3f} a.u. ({how})"
          + ("   [Regenbogen]" if args.rainbow else ""))
    if args.esp_range == "auto":
        print("  ! Fuer den Vergleich mehrerer Molekuele diesen Wert fixieren:")
        print(f"      --esp-range {rng:.3f}")

    # --- Orientierung ---------------------------------------------------
    # Bei mehreren Halogenen blickt die sigma-Ansicht auf das staerkste
    # sigma-Loch - nicht auf das erste Halogen in der Atomliste.
    normal, axis, sigma_axis, center = molecular_frame(
        atoms, halogen_index=loc.get("halogen_index"))  # Bohr (Cube-Einheiten)
    center_ang = center / BOHR_PER_ANGSTROM           # PyMOL rechnet in Angstrom

    views = {
        # Blick senkrecht auf die Ebene; C-X-Achse zeigt nach unten
        "pi":    view_matrix(forward=normal, up=-axis),
        # Blick in der Ebene, senkrecht zur C-X-Achse; C-X-Achse waagerecht
        "edge":  view_matrix(forward=np.cross(normal, axis), up=normal),
        # Blick von aussen entlang der ECHTEN C-X-Achse auf das sigma-Loch.
        # Nicht die in die Ebene projizierte Achse verwenden - bei nicht
        # planaren Molekuelen zeigt die am sigma-Loch vorbei.
        "sigma": view_matrix(forward=sigma_axis, up=normal),
    }
    tilt = float(np.degrees(np.arccos(min(1.0, abs(np.dot(sigma_axis, normal))))))
    tilt = abs(90.0 - tilt)          # Neigung der C-X-Achse gegen die Ebene
    if hals and tilt > 15.0:
        print(f"  Hinweis: die C-X-Achse von {hals[0]['label']} steht "
              f"{tilt:.0f} Grad aus der Ausgleichsebene heraus;")
        print(f"    sigma-Ansicht folgt der echten Bindungsachse, "
              f"pi/edge der Ebene.")
    if args.views:
        views = {k: v for k, v in views.items() if k in args.views}

    # --- PyMOL-Szene ----------------------------------------------------
    cmd.reinitialize()
    cmd.set("auto_zoom", 0)

    cmd.load(args.struct, "mol")
    cmd.load(args.density, "dens")
    cmd.load(args.esp, "esp")

    cmd.hide("everything")
    cmd.show("sticks", "mol")
    cmd.set("stick_radius", 0.10)
    cmd.color("grey20", "mol and elem C")
    cmd.util.cnc("mol")

    cmd.isosurface("surf", "dens", args.iso)
    levels, ramp_colors = ramp_levels(rng, args.rainbow)
    cmd.ramp_new("espramp", "esp", levels, ramp_colors)
    cmd.set("surface_color", "espramp", "surf")
    cmd.disable("espramp")                 # Balken nicht ins Bild rendern

    cmd.set("transparency", args.transparency)
    cmd.set("transparency_mode", 2)
    cmd.set("surface_quality", 1)
    cmd.set("two_sided_lighting", 1)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.15)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)              # keine Perspektive -> vergleichbar

    outdir = args.outdir or "."
    # Eigener Namensanhang, sonst ueberschreibt ein Regenbogenlauf den
    # rot-weiss-blauen Bildersatz desselben Molekuels.
    cmap_tag = "_rainbow" if args.rainbow else ""
    os.makedirs(outdir, exist_ok=True)
    written = []

    for bg in args.backgrounds:
        cmd.bg_color(bg)
        for name, R in views.items():
            v = list(cmd.get_view())
            v[0:9] = [float(x) for x in R.flatten()]
            v[12:15] = [float(x) for x in center_ang]
            cmd.set_view(v)
            # Auf das Molekuel zoomen, NICHT auf "surf": das Isoflaechen-
            # Objekt traegt die Ausdehnung der gesamten Gitterbox mit sich
            # und wuerde das Motiv winzig erscheinen lassen.
            cmd.zoom("mol", args.buffer)

            suffix = f"_{bg}" if len(args.backgrounds) > 1 else ""
            png = os.path.join(outdir, f"{args.prefix}{cmap_tag}_{name}{suffix}.png")
            cmd.ray(args.width, args.height)
            cmd.png(png, dpi=args.dpi)
            written.append(png)
            print(f"    -> {png}")

    # --- Farbskala als eigenes Bild -------------------------------------
    bar = None
    try:
        bar = colorbar(os.path.join(outdir, f"{args.prefix}{cmap_tag}_colorbar.png"),
                       rng, dpi=args.dpi, rainbow=args.rainbow)
        written.append(bar)
        print(f"    -> {bar}")
    except ImportError:
        print("    (matplotlib fehlt - Farbskala wird uebersprungen; "
              "'conda install matplotlib' zum Aktivieren)")

    # --- Protokoll ------------------------------------------------------
    settings = os.path.join(outdir, f"{args.prefix}{cmap_tag}_settings.txt")
    with open(settings, "w", encoding="utf-8") as fh:
        fh.write("Renderparameter (erzeugt von render_esp.py)\n")
        fh.write("=" * 55 + "\n")
        fh.write(f"Struktur          : {args.struct}\n")
        fh.write(f"Dichte-Cube       : {args.density}\n")
        fh.write(f"ESP-Cube          : {args.esp}\n")
        fh.write(f"Gitter            : {dens.shape[0]} x {dens.shape[1]} "
                 f"x {dens.shape[2]}\n")
        fh.write(f"Isowert rho       : {args.iso} a.u.\n")
        fh.write(f"V_S,min           : {vmin:+.5f} a.u. "
                 f"({vmin*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))  auf "
                 f"{loc.get('vmin_atom','?')}\n")
        fh.write(f"V_S,max           : {vmax:+.5f} a.u. "
                 f"({vmax*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))  auf "
                 f"{loc.get('vmax_atom','?')}\n")
        # Eine Zeile pro Halogen, absteigend nach sigma-Loch sortiert.
        for e in hals:
            tag = f"({e['label']})"
            if "sigma_max" not in e:
                fh.write(f"sigma-Loch {tag:<7}: "
                         f"nicht auswertbar (zu wenige Punkte)\n")
                continue
            fh.write(f"sigma-Loch {tag:<7}: "
                     f"{e['sigma_max']:+.5f} a.u. "
                     f"({e['sigma_max']*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))"
                     f"  [{e.get('sigma_method','punktbasiert')}]\n")
            if "belt_min" in e:
                fh.write(f"Guertel    {tag:<7}: "
                         f"{e['belt_min']:+.5f} a.u. "
                         f"({e['belt_min']*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))\n")
        if len(hals) > 1:
            fh.write(f"sigma-Ansicht auf : {hals[0]['label']} "
                     f"(staerkstes sigma-Loch)\n")
        fh.write(f"Gitterabstand     : {spacing:.4f} Bohr\n")
        fh.write(f"Farbskala         : -{rng:.4f} .. +{rng:.4f} a.u. ({how})\n")
        ramp_name = ("Regenbogen (rot-gelb-gruen-cyan-blau)" if args.rainbow
                     else "rot-weiss-blau")
        fh.write(f"Farbrampe         : {ramp_name}\n")
        fh.write(f"Transparenz       : {args.transparency}\n")
        fh.write(f"Hintergrund       : {', '.join(args.backgrounds)}\n")
        fh.write(f"Bildgroesse       : {args.width} x {args.height} px, "
                 f"{args.dpi} dpi\n")
        fh.write(f"Projektion        : orthoskopisch\n")
        fh.write(f"Ansichten         : {', '.join(views.keys())}\n")
    print(f"    -> {settings}")

    return {
        "prefix": args.prefix,
        "struct": args.struct,
        "grid": "x".join(str(v) for v in dens.shape),
        "iso": args.iso,
        "vmin": vmin,
        "vmax": vmax,
        "shell_points": npts,
        "esp_range": rng,
        "esp_range_mode": "auto" if args.esp_range == "auto" else "fixed",
        "colormap": "rainbow" if args.rainbow else "redblue",
        "vmin_atom": loc.get("vmin_atom"),
        "vmax_atom": loc.get("vmax_atom"),
        "halogen": loc.get("halogen"),
        "halogen_atom": loc.get("halogen_atom"),
        "n_halogens": len(hals),
        "sigma_max": loc.get("sigma_max"),
        "sigma_points": loc.get("sigma_points"),
        "sigma_method": loc.get("sigma_method", "punktbasiert"),
        "sigma_angle": loc.get("sigma_angle"),
        "belt_min": loc.get("belt_min"),
        # Alle Halogene, absteigend nach sigma-Loch. Die flachen Felder oben
        # beziehen sich auf halogens[0].
        "halogens": [{k: v for k, v in e.items() if k not in ("pos", "axis")}
                     for e in hals],
        "files": written,
        "settings_file": settings,
    }


def colorbar(path, rng, dpi=300, rainbow=False):
    """Waagerechte Farbskala als separates PNG (braucht matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    cmap = LinearSegmentedColormap.from_list(
        "esp", RAMP_HEX["rainbow" if rainbow else "redblue"])

    # Hoehe grosszuegig plus bbox_inches="tight" beim Speichern: sonst wird die
    # Achsenbeschriftung unten abgeschnitten, was im gerenderten README auffaellt.
    fig = plt.figure(figsize=(4.2, 0.95))
    ax = fig.add_axes([0.06, 0.50, 0.88, 0.26])
    cb = ColorbarBase(ax, cmap=cmap, norm=Normalize(-rng, rng),
                      orientation="horizontal")
    cb.set_label("ESP  /  a.u.", fontsize=9)
    cb.set_ticks([-rng, -rng / 2, 0, rng / 2, rng])
    cb.ax.tick_params(labelsize=8)
    fig.savefig(path, dpi=dpi, transparent=False, facecolor="white",
                bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return path


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def autodetect(pattern_list):
    for pat in pattern_list:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


def main(argv):
    p = argparse.ArgumentParser(
        description="Standardisierte ESP-Bilder aus Cube-Dateien rendern.")
    p.add_argument("--density", default=None, help="Cube der Elektronendichte")
    p.add_argument("--esp", default=None, help="Cube des ESP")
    p.add_argument("--struct", default=None,
                   help="Strukturdatei (.mol/.sdf/.xyz)")
    p.add_argument("--prefix", default=None, help="Praefix der Bildnamen")
    p.add_argument("--outdir", default="images", help="Ausgabeordner")
    p.add_argument("--iso", type=float, default=0.001,
                   help="Isowert der Dichteflaeche in a.u. (Standard 0.001)")
    p.add_argument("--esp-range", default="auto",
                   help="'auto' oder fester Wert in a.u., z.B. 0.03")
    p.add_argument("--transparency", type=float, default=0.15,
                   help="Oberflaechentransparenz 0..1 (Standard 0.15). "
                        "0 = opak, klarste Farben; 0.3+ macht die Profil- "
                        "und Frontalansichten unleserlich, weil man durch "
                        "das ganze Molekuel hindurchschaut.")
    p.add_argument("--backgrounds", nargs="+", default=["white"],
                   help="Hintergrundfarben, z.B. white black")
    p.add_argument("--views", nargs="+", default=None,
                   choices=["pi", "edge", "sigma"],
                   help="Teilmenge der Ansichten (Standard: alle drei)")
    p.add_argument("--width", type=int, default=2000)
    p.add_argument("--height", type=int, default=1600)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--buffer", type=float, default=2.4,
                   help="Rand um das Molekuel in Angstrom")
    p.add_argument("--rainbow", action="store_true",
                   help="Regenbogen-Farbrampe statt rot-weiss-blau. Rot bleibt "
                        "negativ, blau positiv; Gelb/Gruen/Cyan liegen "
                        "dazwischen. Schreibt einen eigenen Bildersatz "
                        "<prefix>_rainbow_*.png")
    p.add_argument("--no-color", action="store_true",
                   help="plain terminal output without ANSI colours (same effect as "
                        "setting the NO_COLOR environment variable)")
    args = p.parse_args(argv)

    if args.no_color:
        ansi.disable()

    args.density = args.density or autodetect(["td.cube", "*dens*.cube"])
    args.esp = args.esp or autodetect(["tp.cube", "*esp*.cube", "*pot*.cube"])
    args.struct = args.struct or autodetect(
        ["*.mol", "*.sdf", "*.xyz"])

    missing = [n for n, v in (("--density", args.density),
                              ("--esp", args.esp),
                              ("--struct", args.struct)) if not v]
    if missing:
        raise SystemExit("Fehlende Eingaben: " + ", ".join(missing))

    if not args.prefix:
        base = os.path.splitext(os.path.basename(args.struct))[0]
        args.prefix = base.split("_")[0] or "molecule"

    print("=" * 70)
    print("render_esp.py - standardisierte ESP-Bilder")
    print("=" * 70)
    print(f"  Struktur : {args.struct}")
    print(f"  Dichte   : {args.density}")
    print(f"  ESP      : {args.esp}")
    print(f"  Praefix  : {args.prefix}")

    render_all(args)
    print("Fertig.")
    return 0


# Ausfuehren, sobald das Skript NICHT als Modul importiert wird.
#
# Warum nicht das uebliche  if __name__ == "__main__"  ?
# PyMOL fuehrt uebergebene .py-Dateien mit exec() in einem eigenen Namensraum
# aus, in dem __name__ eben nicht "__main__" ist. Mit der Standardabfrage
# passiert bei  pymol -cq render_esp.py -- ...  schlicht gar nichts:
# das Skript wird gelesen, alle Funktionen werden definiert, und dann ist
# Schluss - ohne Fehlermeldung. Genau diese stille Nicht-Ausfuehrung ist
# schwer zu diagnostizieren, deshalb hier die umgekehrte Abfrage.
if __name__ != "render_esp":
    _argv = sys.argv[1:]
    if "--" in _argv:                  # Aufruf ueber: pymol -cq skript.py -- ...
        _argv = _argv[_argv.index("--") + 1:]
    else:
        # PyMOL schiebt beim Start eigene Argumente in sys.argv. Alles vor
        # der Skriptdatei wegwerfen, damit argparse nicht darueber stolpert.
        for _i, _a in enumerate(_argv):
            if _a.endswith("render_esp.py"):
                _argv = _argv[_i + 1:]
                break
    main(_argv)
