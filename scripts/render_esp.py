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
``tp.cube`` und einer Struktur (.mol/.sdf/.pdb/.xyz).


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
import math
import os
import sys

import numpy as np

import ansi


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


def esp_statistics(density, esp, iso=0.001, tol_factor=0.12):
    """ESP-Kennzahlen auf der rho=iso-Schale.

    Liefert (V_min, V_max, anzahl_punkte) in atomaren Einheiten.
    ``tol_factor`` legt die Schalendicke relativ zum Isowert fest.
    """
    mask = np.abs(density - iso) < iso * tol_factor
    if mask.sum() < 50:                        # Schale zu duenn -> aufweiten
        mask = np.abs(density - iso) < iso * 0.30
    if mask.sum() == 0:
        return None, None, 0
    shell = esp[mask]
    return float(shell.min()), float(shell.max()), int(mask.sum())


def shell_points(density, esp, origin, voxel, iso=0.001, tol_factor=0.12):
    """Koordinaten und ESP-Werte der Gitterpunkte auf der rho=iso-Schale.

    Nur die Schalenpunkte werden materialisiert, nicht das ganze Gitter -
    bei 251^3 waere ein volles Koordinatenfeld sonst mehrere hundert MB.
    """
    mask = np.abs(density - iso) < iso * tol_factor
    if mask.sum() < 50:
        mask = np.abs(density - iso) < iso * 0.30
    idx = np.argwhere(mask)                       # (N, 3) Gitterindizes
    pos = origin + idx @ voxel                    # (N, 3) kartesisch, Bohr
    return pos, esp[mask]


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
      hydro   Umgebung der Wasserstoffatome

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

    # --- Halogenregionen ------------------------------------------------
    hal = [i for i, z in enumerate(znums) if z in HALOGENS]
    if not hal:
        return out

    hi = hal[0]
    carbons = [i for i, z in enumerate(znums) if z == 6]
    if not carbons:
        return out
    ci = carbons[int(np.argmin(np.linalg.norm(coords[carbons] - coords[hi],
                                              axis=1)))]
    axis = coords[hi] - coords[ci]
    axis = axis / np.linalg.norm(axis)             # C -> X, zeigt zum sigma-Loch

    rel = pos - coords[hi]
    r = np.linalg.norm(rel, axis=1)
    r[r == 0] = 1e-9
    cos = (rel @ axis) / r

    cap = cos > cone_cos
    if cap.sum() >= 5:
        out["sigma_max"] = float(vals[cap].max())
        out["sigma_points"] = int(cap.sum())
        # Wie weit sitzt der gefundene Maximalpunkt von der C-X-Achse entfernt?
        # Das sigma-Loch ist ein Gipfel AUF der Achse; wird die vom Gitter nicht
        # getroffen, misst man die Flanke und unterschaetzt den Wert. Die reine
        # Punktzahl in der Kappe reicht als Kriterium nicht aus: bei Brombenzol
        # lagen 144 Punkte in der Kappe, der beste davon aber 1.14 Bohr neben
        # der Achse - Ergebnis +0.0126 statt +0.0175 auf dem feinen Gitter.
        j = int(np.argmax(np.where(cap, vals, -np.inf)))
        out["sigma_offaxis"] = float(r[j] * np.sqrt(max(0.0, 1 - cos[j] ** 2)))
        # Die sigma-Kappe ist ein kleiner Ausschnitt der Oberflaeche. Auf einem
        # groben Gitter liegen dort nur wenige Punkte, und das lokale Maximum
        # wird dann systematisch unterschaetzt.
        out["sigma_sparse"] = bool(cap.sum() < 30)
        r_cap = float(r[cap].mean())
        belt = (np.abs(cos) < belt_cos) & (r < belt_factor * r_cap)
        if belt.sum() >= 5:
            out["belt_min"] = float(vals[belt].min())
            out["belt_points"] = int(belt.sum())

    out["halogen"] = HALOGENS[int(znums[hi])]
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
    """Gleichmaessig verteilte Richtungen in einer Kappe um ``axis``.

    Fibonacci-Spirale auf der Kugelkappe - liefert eine gleichmaessige Belegung
    ohne Haeufung an der Achse, wie sie bei Kugelkoordinaten auftraete.
    """
    k = np.arange(n) + 0.5
    cosv = 1.0 - (1.0 - cone_cos) * k / n          # cone_cos .. 1
    phi = np.pi * (1 + 5 ** 0.5) * k
    sinv = np.sqrt(np.maximum(0.0, 1 - cosv ** 2))

    # orthonormale Basis um axis
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, axis)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    return (cosv[:, None] * axis
            + (sinv * np.cos(phi))[:, None] * e1
            + (sinv * np.sin(phi))[:, None] * e2)


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

    Rueckgabe: dict mit sigma_max, belt_min und dem Winkel des Maximums zur
    Achse; None, wenn kein Halogen vorhanden ist.
    """
    diag = np.diag(voxel)
    if not np.allclose(voxel, np.diag(diag)):
        return None                       # nicht achsparallel - Fallback
    delta = diag

    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])
    hal = [i for i, z in enumerate(znums) if z in HALOGENS]
    carbons = [i for i, z in enumerate(znums) if z == 6]
    if not hal or not carbons:
        return None

    hi = hal[0]
    ci = carbons[int(np.argmin(np.linalg.norm(coords[carbons] - coords[hi],
                                              axis=1)))]
    axis = coords[hi] - coords[ci]
    axis = axis / np.linalg.norm(axis)

    dirs = _cone_directions(axis, cone_cos, n_rays)
    radii = np.arange(1.0, r_max, dr)

    best_v, best_cos = None, None
    for d in dirs:
        pts = coords[hi] + radii[:, None] * d[None, :]
        rho = _trilinear(density, origin, delta, pts)
        # aeusserster Schnittpunkt: von aussen nach innen der erste Wert >= iso
        idx = np.nonzero(rho >= iso)[0]
        if idx.size == 0:
            continue
        j = idx[-1]
        if j + 1 >= len(radii):
            continue
        # lineare Interpolation des Radius am Isowert
        r0, r1 = radii[j], radii[j + 1]
        y0, y1 = rho[j], rho[j + 1]
        rs = r0 + (iso - y0) * (r1 - r0) / (y1 - y0) if y1 != y0 else r0
        v = float(_trilinear(esp, origin, delta,
                             (coords[hi] + rs * d)[None, :])[0])
        if best_v is None or v > best_v:
            best_v, best_cos = v, float(np.dot(d, axis))

    if best_v is None:
        return None
    return {"sigma_max": best_v,
            "sigma_angle": float(np.degrees(np.arccos(min(1.0, best_cos)))),
            "sigma_method": "interpoliert"}


def nice_range(vmin, vmax, step=0.005):
    """Symmetrischer, auf ``step`` aufgerundeter Bereich."""
    amp = max(abs(vmin), abs(vmax))
    return math.ceil(amp / step) * step


# ----------------------------------------------------------------------------
# Orientierung aus der Geometrie
# ----------------------------------------------------------------------------

HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}
BOHR_PER_ANGSTROM = 1.8897259886

Z_SYMBOL = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na",
            12: "Mg", 14: "Si", 15: "P", 16: "S", 17: "Cl", 19: "K",
            20: "Ca", 26: "Fe", 29: "Cu", 30: "Zn", 34: "Se", 35: "Br",
            53: "I"}


def z_symbol(z):
    return Z_SYMBOL.get(int(z), f"Z{int(z)}")


def molecular_frame(atoms):
    """Bestimmt ein reproduzierbares Molekuelkoordinatensystem.

    Rueckgabe: (normal, axis, center)
      normal  Flaechennormale (kleinste Traegheitsausdehnung, nur Schweratome)
      axis    C->Halogen-Achse; falls kein Halogen: laengste Hauptachse
      center  geometrischer Mittelpunkt aller Atome
    Alle Vektoren normiert, Koordinaten in denselben Einheiten wie ``atoms``.
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
        hi = hal_idx[0]
        carbons = [i for i, z in enumerate(znums) if z == 6]
        if carbons:
            d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
            ci = carbons[int(np.argmin(d))]
            axis = coords[hi] - coords[ci]           # C -> X, zeigt zum sigma-Loch

    if axis is None:
        axis = long_axis.copy()

    axis = axis / np.linalg.norm(axis)
    normal = normal / np.linalg.norm(normal)

    # axis exakt senkrecht zur Normalen machen (numerisches Aufraeumen)
    axis = axis - normal * float(np.dot(axis, normal))
    if np.linalg.norm(axis) < 1e-6:
        axis = long_axis
    axis = axis / np.linalg.norm(axis)

    return normal, axis, center


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
    if ray:
        loc.update(ray)
    spacing = float(np.max(np.abs(np.diag(voxel))))

    def _fmt(v):
        return (f"{v:+.4f} a.u.  = {v*2625.4996:+7.1f} kJ/(mol*e)"
                f"  = {v*627.5095:+6.1f} kcal/(mol*e)")

    print(f"  ESP auf der rho={args.iso}-Schale ({npts} Punkte):")
    print(f"    V_S,min = {_fmt(vmin)}   auf "
          f"{ansi.atom_label(loc.get('vmin_atom', '?'))}")
    print(f"    V_S,max = {_fmt(vmax)}   auf "
          f"{ansi.atom_label(loc.get('vmax_atom', '?'))}")
    if "sigma_max" in loc:
        print(f"  Lokal am Halogen "
              f"({ansi.element(loc.get('halogen', '?'))}):")
        tag = loc.get("sigma_method", "punktbasiert")
        extra = (f"   [{tag}, {loc['sigma_angle']:.1f} Grad zur Achse]"
                 if "sigma_angle" in loc
                 else f"   [{tag}, {loc.get('sigma_points', 0)} Punkte]")
        print(f"    sigma-Loch  = {_fmt(loc['sigma_max'])}{extra}")
        # Auch das Strahlverfahren kann die Isoflaeche auf einem groben Gitter
        # nur so genau lokalisieren, wie die Dichte dort aufgeloest ist.
        if spacing > 0.30:
            print(f"    ! Gitterabstand {spacing:.2f} Bohr - fuer einen "
                  f"belastbaren sigma-Loch-Wert zu grob;")
            print(f"      erwartungsgemaess einige Prozent zu niedrig. "
                  f"Feiner rechnen (kleineres --stride).")
        if "belt_min" in loc:
            print(f"    Guertel     = {_fmt(loc['belt_min'])}"
                  f"   [{loc['belt_points']} Punkte]")
        # Hinweis, dass V_S,max nicht auf dem Halogen liegt: bewusst nicht
        # ausgegeben. Bei Arylhalogeniden trifft das praktisch immer zu, die
        # Meldung waere also bei jedem Molekuel identisch und damit wertlos.
        # Die Information steckt bereits in der Ortsangabe hinter V_S,max
        # ("auf H5") und im separat ausgewiesenen sigma-Loch. Erklaerung dazu
        # in der README, Abschnitt "Which number describes the sigma-hole".
    print(f"  Farbskala: +/- {rng:.3f} a.u. ({how})")
    if args.esp_range == "auto":
        print("  ! Fuer den Vergleich mehrerer Molekuele diesen Wert fixieren:")
        print(f"      --esp-range {rng:.3f}")

    # --- Orientierung ---------------------------------------------------
    normal, axis, center = molecular_frame(atoms)     # in Bohr (Cube-Einheiten)
    center_ang = center / BOHR_PER_ANGSTROM           # PyMOL rechnet in Angstrom

    views = {
        # Blick senkrecht auf die Ebene; C-X-Achse zeigt nach unten
        "pi":    view_matrix(forward=normal, up=-axis),
        # Blick in der Ebene, senkrecht zur C-X-Achse; C-X-Achse waagerecht
        "edge":  view_matrix(forward=np.cross(normal, axis), up=normal),
        # Blick von aussen entlang der C-X-Achse auf das sigma-Loch
        "sigma": view_matrix(forward=axis, up=normal),
    }
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
    cmd.ramp_new("espramp", "esp", [-rng, 0.0, rng], ["red", "white", "blue"])
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
            png = os.path.join(outdir, f"{args.prefix}_{name}{suffix}.png")
            cmd.ray(args.width, args.height)
            cmd.png(png, dpi=args.dpi)
            written.append(png)
            print(f"    -> {png}")

    # --- Farbskala als eigenes Bild -------------------------------------
    bar = None
    try:
        bar = colorbar(os.path.join(outdir, f"{args.prefix}_colorbar.png"),
                       rng, dpi=args.dpi)
        written.append(bar)
        print(f"    -> {bar}")
    except ImportError:
        print("    (matplotlib fehlt - Farbskala wird uebersprungen; "
              "'conda install matplotlib' zum Aktivieren)")

    # --- Protokoll ------------------------------------------------------
    settings = os.path.join(outdir, f"{args.prefix}_settings.txt")
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
                 f"({vmin*627.5095:+.2f} kcal/(mol*e))  auf "
                 f"{loc.get('vmin_atom','?')}\n")
        fh.write(f"V_S,max           : {vmax:+.5f} a.u. "
                 f"({vmax*627.5095:+.2f} kcal/(mol*e))  auf "
                 f"{loc.get('vmax_atom','?')}\n")
        if "sigma_max" in loc:
            fh.write(f"sigma-Loch ({loc.get('halogen','?'):<2})   : "
                     f"{loc['sigma_max']:+.5f} a.u. "
                     f"({loc['sigma_max']*627.5095:+.2f} kcal/(mol*e))"
                     f"  [{loc.get('sigma_method','punktbasiert')}]\n")
        fh.write(f"Gitterabstand     : {spacing:.4f} Bohr\n")
        if "belt_min" in loc:
            fh.write(f"Halogenguertel    : {loc['belt_min']:+.5f} a.u. "
                     f"({loc['belt_min']*627.5095:+.2f} kcal/(mol*e))\n")
        fh.write(f"Farbskala         : -{rng:.4f} .. +{rng:.4f} a.u. ({how})\n")
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
        "vmin_atom": loc.get("vmin_atom"),
        "vmax_atom": loc.get("vmax_atom"),
        "halogen": loc.get("halogen"),
        "sigma_max": loc.get("sigma_max"),
        "sigma_points": loc.get("sigma_points"),
        "sigma_method": loc.get("sigma_method", "punktbasiert"),
        "sigma_angle": loc.get("sigma_angle"),
        "belt_min": loc.get("belt_min"),
        "files": written,
        "settings_file": settings,
    }


def colorbar(path, rng, dpi=300):
    """Waagerechte Farbskala als separates PNG (braucht matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    cmap = LinearSegmentedColormap.from_list(
        "esp", ["#d40000", "#ffffff", "#0030d4"])

    fig = plt.figure(figsize=(4.2, 0.75))
    ax = fig.add_axes([0.06, 0.45, 0.88, 0.30])
    cb = ColorbarBase(ax, cmap=cmap, norm=Normalize(-rng, rng),
                      orientation="horizontal")
    cb.set_label("ESP  /  a.u.", fontsize=9)
    cb.set_ticks([-rng, -rng / 2, 0, rng / 2, rng])
    cb.ax.tick_params(labelsize=8)
    fig.savefig(path, dpi=dpi, transparent=False, facecolor="white")
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
                   help="Strukturdatei (.mol/.sdf/.pdb/.xyz)")
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
    p.add_argument("--no-color", action="store_true",
                   help="plain output without ANSI colours (same effect as "
                        "setting the NO_COLOR environment variable)")
    args = p.parse_args(argv)

    if args.no_color:
        ansi.disable()

    args.density = args.density or autodetect(["td.cube", "*dens*.cube"])
    args.esp = args.esp or autodetect(["tp.cube", "*esp*.cube", "*pot*.cube"])
    args.struct = args.struct or autodetect(
        ["*.mol", "*.sdf", "*.pdb", "*.xyz"])

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
