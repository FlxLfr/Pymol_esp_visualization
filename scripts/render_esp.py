#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_esp.py
=============

Produces a standardised set of ESP images from a pair of cube files (electron
density + electrostatic potential), fully automatic and without a single mouse
click.

The ESP is mapped onto the isosurface of the electron density at
rho = 0.001 a.u. (the Politzer/Murray convention) and rendered from three
fixed viewing directions:

    pi      perpendicular to the molecular plane -> shows the pi system
    edge    in the molecular plane               -> profile, C-X axis horizontal
    sigma   along the C-X axis from outside      -> shows the sigma hole head on

The orientation is computed from the geometry (inertial axes + the
carbon-halogen axis), NOT via PyMOL's ``orient``. That way different molecules
reproducibly come out in the same alignment.


Call
----
    pymol -cq render_esp.py -- --density td.cube --esp tp.cube \\
                               --struct brombenzol_aro_opti.mol \\
                               --prefix brombenzol

Without arguments the script looks in the current folder for ``td.cube``,
``tp.cube`` and a structure (.mol/.sdf/.xyz).


Colour scale
------------
By default the ESP range is determined from the data *on the isosurface* and
rounded up symmetrically to a round value. The value used is written to the log
and to the file ``<prefix>_settings.txt``.

!! For a direct comparison of several molecules the scale has to be fixed.
   Run all molecules once with --esp-range auto, note the largest value
   reported, and then render them all again with --esp-range <value>.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

import ansi
import xyzToCube                    # element list (see Z_SYMBOL) and the shell
# The shell evaluation (rho = iso) and the colour scale live in xyzToCube.py,
# because the converter needs them for its own esp.pml as well. Both scripts
# therefore use the same definition of "on the isosurface".
from xyzToCube import esp_statistics, nice_range, shell_mask, SHELL_TOL_FACTOR
from constants import BOHR_PER_ANGSTROM, HARTREE_TO_KCAL, HARTREE_TO_KJ


# ----------------------------------------------------------------------------
# Reading the cube (for the statistics only; PyMOL loads the files itself)
# ----------------------------------------------------------------------------

def read_cube(path):
    """Reads a Gaussian cube. Returns (values3d, atoms, origin, voxel)."""
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
        raise ValueError(f"{path}: expected {int(np.prod(n))} values, "
                         f"read {values.size}")
    return values.reshape(n), atoms, origin, voxel


def shell_points(density, esp, origin, voxel, iso=0.001,
                 tol_factor=SHELL_TOL_FACTOR):
    """Coordinates and ESP values of the grid points on the rho = iso shell.

    Only the shell points are materialised, not the whole grid - at 251^3 a
    full coordinate array would otherwise be several hundred MB.
    """
    mask = shell_mask(density, iso, tol_factor)
    idx = np.argwhere(mask)                       # (N, 3) grid indices
    pos = origin + idx @ voxel                    # (N, 3) cartesian, Bohr
    return pos, esp[mask]


def halogen_axes(atoms):
    """All halogens with their C-X axis.

    A molecule can carry more than one halogen - each is evaluated separately.

    Returns a list of dicts with
      index   0-based atom index of the halogen
      symbol  element symbol
      label   symbol + 1-based number, e.g. "Cl12" (as in the output)
      pos     coordinates of the halogen
      axis    normalised C->X axis (points at the sigma hole)
      r_limit the largest distance at which the halogen's own surface can lie
              (Bohr) - see below
    Halogens without a bonded carbon within reach are skipped.

    Why ``r_limit``: in folded molecules the cone around the C-X axis does not
    point into empty space but at another part of the molecule. Without a
    distance limit that part's surface is measured as well. Triazolam is
    exactly such a case - the cone around Cl21 hits the methyl group on the
    triazole ring, and that returned a "sigma hole" of +18.8 instead of +10.5
    kcal/(mol*e). The rho = 0.001 surface of a halogen sits at about 1.1 to 1.2
    vdW radii; the factor 1.6 leaves ample room and excludes everything beyond.
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
    """Region-resolved surface extrema.

    Why this is needed: the *global* V_S,max of an aryl halide surface sits
    almost always on the ring hydrogens, not on the halogen. For bromobenzene
    and iodobenzene the global value therefore returns the same C-H twice
    (+0.031 a.u.), while the sigma holes, which are the actual subject, differ
    by almost a factor of two. Any comparison of halogen-bond donors needs the
    *local* maximum on the halogen.

    Regions:
      sigma   cap around the extended C-X axis (opening angle from cone_cos)
      belt    belt perpendicular to it, close to the halogen

    Returns a dict with the values in atomic units; the halogen-related entries
    are absent if the molecule contains no halogen.
    """
    out = {}
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])

    # --- name the region of the global extremum ------------------------
    for tag, i in (("vmax", int(np.argmax(vals))), ("vmin", int(np.argmin(vals)))):
        d = np.linalg.norm(coords - pos[i], axis=1)
        j = int(np.argmin(d))
        out[f"{tag}_atom"] = f"{z_symbol(int(znums[j]))}{j + 1}"

    # --- halogen regions, one after the other ----------------------------
    entries = []
    for hal in halogen_axes(atoms):
        e = {k: hal[k] for k in ("index", "symbol", "label")}
        axis = hal["axis"]

        rel = pos - hal["pos"]
        r = np.linalg.norm(rel, axis=1)
        r[r == 0] = 1e-9
        cos = (rel @ axis) / r

        # Distance limit: otherwise, in folded molecules, points on a
        # completely different part of the molecule count towards the cap (see
        # halogen_axes).
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
    """Sort the halogen entries by sigma hole, descending.

    Entries without an evaluable sigma hole go to the end. The order determines
    the output and - via the first entry - the orientation of the sigma view.
    """
    return sorted(entries,
                  key=lambda e: -e.get("sigma_max", -np.inf))


def promote_primary(out):
    """Also place the values of the strongest sigma hole flat into ``out``.

    That keeps callers working which expect only ONE value (the CSV column
    ``sigma_hole_au``, the summary table); with a single halogen the result is
    bit-identical to before.
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
    """Trilinear interpolation on an axis-aligned, regular grid."""
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
    """Directions in a cap around ``axis``; the first one is the axis itself.

    The rest is a Fibonacci spiral on the spherical cap - even coverage without
    the crowding near the axis that spherical coordinates would produce.

    Why the axis is added separately: the ``+ 0.5`` in ``k`` is the midpoint
    rule for an EQUAL-AREA distribution - every ray sits in the middle of a
    ring of the same size. That is right when averaging over the cap. But we
    are looking for a MAXIMUM, and for the sigma hole that sits exactly on the
    axis. With the offset the innermost ray was 1.281 degrees away from it, the
    axis itself was never evaluated, and every axially symmetric molecule
    stubbornly reported "1.3 degrees" - a lower bound of the sampling grid, not
    a measurement. The error in the value was small (4-bromoacetophenone: 0.008
    kcal/(mol*e)), but the report was misleading.
    """
    axis = np.asarray(axis, dtype=float)

    # n-1 spiral directions; the axis is prepended
    m = max(1, n - 1)
    k = np.arange(m) + 0.5
    cosv = 1.0 - (1.0 - cone_cos) * k / m          # cone_cos .. 1
    phi = np.pi * (1 + 5 ** 0.5) * k
    sinv = np.sqrt(np.maximum(0.0, 1 - cosv ** 2))

    # orthonormal basis around axis
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
    """The sigma hole, independent of how the grid is aligned.

    The problem with the point-based evaluation: the sigma hole is a peak *on*
    the C-X axis. Whether a grid point happens to sit there AND inside the thin
    rho = iso shell at the same time is luck. For bromobenzene the best point
    on the 126^3 grid lay 1.14 Bohr off the axis - the value came out 28 % too
    low, although 144 points lay inside the cap.

    Instead, here: rays from the halogen into a cap around the axis, on every
    ray the radius at which rho crosses the isosurface, and V evaluated there -
    both trilinearly interpolated. The result no longer depends on where the
    grid points happen to lie.

    If the molecule carries several halogens, each is sampled separately.

    Returns a dict {atom index: {sigma_max, sigma_angle, sigma_method}}; an
    empty dict if no halogen can be evaluated.
    """
    diag = np.diag(voxel)
    if not np.allclose(voxel, np.diag(diag)):
        return {}                         # not axis-aligned - fall back
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
            # INNERMOST crossing, from inside outwards: the ray starts deep
            # inside the halogen's density, and the first drop below iso is its
            # own surface. The outermost crossing used to be taken; in folded
            # molecules the ray then dives into another part of the molecule
            # behind it and measures that surface instead.
            if rho[0] < iso:
                continue                  # the ray already starts outside
            below = np.nonzero(rho < iso)[0]
            if below.size == 0:
                continue                  # surface not hit within r_limit
            j = below[0] - 1
            if j < 0 or j + 1 >= len(radii):
                continue
            # linear interpolation of the radius at the isovalue
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
            "sigma_method": "interpolated"}

    return result


# ----------------------------------------------------------------------------
# Orientation from the geometry
# ----------------------------------------------------------------------------

HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}

# ----------------------------------------------------------------------------
# Colour ramps
#
# Both keep the same convention: RED = negative, BLUE = positive, and the
# middle of the scale is V = 0. The rainbow ramp only pushes yellow/green/cyan
# in between instead of passing through white. That keeps the images of both
# ramps comparable at the ends - were the rainbow scale reversed (blue =
# negative, as in some programs), the two sets could not be laid side by side.
#
# What the rainbow ramp is good for: red-white-blue has almost no colour
# resolution in the middle, and weakly polar regions all look equally white.
# The rainbow resolves exactly there. The price: it is not perceptually uniform
# and creates edges where there are none - for a quantitative statement
# red-white-blue remains the more honest depiction.
# ----------------------------------------------------------------------------

# Margin around the molecule in the image, in Angstrom.
#
# The zoom is on the ATOMS, not on the isosurface: the surface object carries
# the extent of the whole grid box with it and would leave the molecule looking
# like a postage stamp in the middle of the image. The room the surface needs
# beyond the nuclei therefore has to be added here - the rho = 0.001 surface
# sits about 1.7 to 2.1 Angstrom outside the outermost nuclei, and 2.4 covers
# that with a narrow margin.
#
# This used to be an option (--buffer). It was taken out because no measured
# value depends on it, and a different margin only produces an image set that
# can no longer be laid beside the others.
BUFFER_ANGSTROM = 2.4

RAMP_PYMOL = {
    "redblue": ["red", "white", "blue"],
    "rainbow": ["red", "yellow", "green", "cyan", "blue"],
}
RAMP_HEX = {
    "redblue": ["#d40000", "#ffffff", "#0030d4"],
    "rainbow": ["#d40000", "#f0e000", "#00a000", "#00c8d4", "#0030d4"],
}


def ramp_levels(rng, rainbow=False):
    """Anchor levels and colours for ``cmd.ramp_new``.

    Returns (levels, colors) - equal length, symmetric about 0.
    """
    name = "rainbow" if rainbow else "redblue"
    colors = RAMP_PYMOL[name]
    n = len(colors)
    levels = [-rng + 2.0 * rng * i / (n - 1) for i in range(n)]
    return levels, colors

# van der Waals radii after Bondi (J. Phys. Chem. 1964, 68, 441), Angstrom.
# Only as an order of magnitude for the distance limit of the sigma-hole
# search.
VDW_ANGSTROM = {9: 1.47, 17: 1.75, 35: 1.85, 53: 1.98}

# Reverse lookup, atomic number -> symbol, derived from the same list with
# which xyzToCube.py translates in the other direction. There used to be a
# second, hand-maintained table of 20 entries here: a molecule with an element
# missing from it converted cleanly and was then labelled "Z13" instead of
# "Al13". One source, both directions.
Z_SYMBOL = {i + 1: sym for i, sym in enumerate(xyzToCube.ELEMENTS)}


def z_symbol(z):
    return Z_SYMBOL.get(int(z), f"Z{int(z)}")


def molecular_frame(atoms, halogen_index=None):
    """Determines a reproducible molecular frame.

    With several halogens, ``halogen_index`` selects which one fixes the axis -
    and therefore which sigma hole the sigma view looks at. Without it the
    first halogen in the atom list is taken; render_all passes the halogen with
    the strongest sigma hole.

    Returns (normal, axis, sigma_axis, center)
      normal      surface normal (smallest inertial extent, heavy atoms)
      axis        the C->halogen axis PROJECTED INTO THE PLANE. Together with
                  ``normal`` it spans a clean right-handed frame and orients
                  the pi and edge views.
      sigma_axis  the TRUE C->halogen axis, unprojected
      center      geometric centre of all atoms
    All vectors normalised, coordinates in the same units as ``atoms``.

    Why two axes: for planar molecules both are identical and the projection is
    pure rounding cosmetics there. For non-planar molecules it really does turn
    the axis away - in triazolam the C-Cl21 bond points 42.9 degrees out of the
    best-fit plane, and the sigma view accordingly looked past the sigma hole
    by 42.9 degrees. For the sigma view, therefore, use ``sigma_axis``.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    znums = np.array([a[0] for a in atoms])
    center = coords.mean(axis=0)

    heavy = coords[znums > 1]
    if len(heavy) < 3:
        heavy = coords
    centered = heavy - heavy.mean(axis=0)

    # principal axes via singular value decomposition
    _, sing, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[2]                                   # smallest extent
    long_axis = vt[0]                                # largest extent

    # find the C-halogen axis
    axis = None
    hal_idx = [i for i, z in enumerate(znums) if z in HALOGENS]
    if hal_idx:
        hi = halogen_index if halogen_index in hal_idx else hal_idx[0]
        carbons = [i for i, z in enumerate(znums) if z == 6]
        if carbons:
            d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
            ci = carbons[int(np.argmin(d))]
            axis = coords[hi] - coords[ci]           # C -> X, points at the sigma hole

    if axis is None:
        axis = long_axis.copy()

    axis = axis / np.linalg.norm(axis)
    normal = normal / np.linalg.norm(normal)
    sigma_axis = axis.copy()                 # unchanged, for the sigma view

    # project axis into the plane - only for pi and edge, which need a
    # right-handed frame together with the normal.
    axis = axis - normal * float(np.dot(axis, normal))
    if np.linalg.norm(axis) < 1e-6:
        axis = long_axis
    axis = axis / np.linalg.norm(axis)

    return normal, axis, sigma_axis, center


def view_matrix(forward, up):
    """Rotation matrix for PyMOL's set_view.

    ``forward`` points from the molecule to the camera, ``up`` points up in the
    image. The rows of the matrix are the camera basis vectors in world space.
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
    # In set_view PyMOL expects the matrix whose COLUMNS are the camera basis
    # vectors in world space - that is, the transpose of the row form. Checked
    # empirically (see the SOP, section on views).
    return np.array([right, up, z]).T


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------

def ensure_pymol():
    """Returns ``cmd``; starts PyMOL headless if it is not running yet.

    That way the script runs both as ``python render_esp.py ...`` and as
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

    # --- data for the statistics ----------------------------------------
    dens, atoms, origin, voxel = read_cube(args.density)
    esp, _, _, _ = read_cube(args.esp)
    if dens.shape != esp.shape:
        raise SystemExit("Density and ESP cube are on different grids.")

    vmin, vmax, npts = esp_statistics(dens, esp, iso=args.iso)
    if npts == 0:
        raise SystemExit(f"No grid points found at rho = {args.iso}. "
                         f"Check the isovalue.")

    if args.esp_range == "auto":
        rng = nice_range(vmin, vmax)
        how = "automatic, from the data"
    else:
        rng = float(args.esp_range)
        how = "fixed"

    # Region-resolved values: for aryl halides the global maximum sits on the
    # ring hydrogens, not on the halogen.
    pos, vals = shell_points(dens, esp, origin, voxel, iso=args.iso)
    loc = local_extrema(pos, vals, atoms)

    # The sigma hole is preferably determined by rays with interpolation. The
    # point-based variant depends on whether a grid point happens to lie near
    # the C-X axis AND inside the thin isosurface shell; for bromobenzene it
    # gives +7.9 instead of +10.1 kcal/(mol*e) on the same grid.
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

    print(f"  ESP on the rho = {args.iso} shell ({npts} points):")
    print(f"    V_S,min = {_fmt(vmin)}   on "
          f"{ansi.atom_label(loc.get('vmin_atom', '?'))}")
    print(f"    V_S,max = {_fmt(vmax)}   on "
          f"{ansi.atom_label(loc.get('vmax_atom', '?'))}")
    # One block per halogen. With exactly one halogen the output is the same
    # as before; from two on, each gets its own line, and the halogen the sigma
    # view looks at is marked.
    for e in hals:
        if "sigma_max" not in e:
            print(f"  Local at the halogen ({ansi.element(e['symbol'])}"
                  f" {e['index'] + 1}): no evaluable sigma hole "
                  f"(too few surface points in the cap)")
            continue
        head = (f"  Local at the halogen ({ansi.element(e['symbol'])}):"
                if len(hals) == 1
                else f"  Local at {ansi.atom_label(e['label'])}:")
        if len(hals) > 1 and e is hals[0]:
            head += "   <- axis of the sigma view"
        print(head)
        tag = e.get("sigma_method", "point-based")
        extra = (f"   [{tag}, {e['sigma_angle']:.1f} degrees off the axis]"
                 if "sigma_angle" in e
                 else f"   [{tag}, {e.get('sigma_points', 0)} points]")
        print(f"    sigma hole  = {_fmt(e['sigma_max'])}{extra}")
        if "belt_min" in e:
            print(f"    belt        = {_fmt(e['belt_min'])}"
                  f"   [{e['belt_points']} points]")
    if hals and any("sigma_max" in e for e in hals):
        # The ray method, too, can only locate the isosurface on a coarse grid
        # as precisely as the density is resolved there.
        if spacing > 0.30:
            print(f"    ! grid spacing {spacing:.2f} Bohr - too coarse for a "
                  f"trustworthy sigma-hole value;")
            print(f"      expect it to be a few per cent low. Compute finer "
                  f"(a smaller --stride).")
        # A note that V_S,max does not sit on the halogen is deliberately not
        # printed. For aryl halides that is practically always the case, so the
        # message would be identical for every molecule and therefore
        # worthless. The information is already in the location behind V_S,max
        # ("on H5") and in the separately reported sigma hole. Explained in
        # docs/ESP_Visualization_Background.docx, section 2.1
        # "Which number describes the sigma-hole - Not V_S,max".
    print(f"  Colour scale: +/- {rng:.3f} a.u. ({how})"
          + ("   [rainbow]" if args.rainbow else ""))
    if args.esp_range == "auto":
        print("  ! To compare several molecules, fix this value:")
        print(f"      --esp-range {rng:.3f}")

    # --- orientation ----------------------------------------------------
    # With several halogens the sigma view looks at the strongest sigma hole -
    # not at the first halogen in the atom list.
    normal, axis, sigma_axis, center = molecular_frame(
        atoms, halogen_index=loc.get("halogen_index"))  # Bohr (cube units)
    center_ang = center / BOHR_PER_ANGSTROM           # PyMOL works in Angstrom

    views = {
        # looking perpendicular onto the plane; C-X axis points down
        "pi":    view_matrix(forward=normal, up=-axis),
        # looking in the plane, perpendicular to the C-X axis; C-X horizontal
        "edge":  view_matrix(forward=np.cross(normal, axis), up=normal),
        # looking in from outside along the TRUE C-X axis onto the sigma hole.
        # Do not use the axis projected into the plane - for non-planar
        # molecules it points past the sigma hole.
        "sigma": view_matrix(forward=sigma_axis, up=normal),
    }
    tilt = float(np.degrees(np.arccos(min(1.0, abs(np.dot(sigma_axis, normal))))))
    tilt = abs(90.0 - tilt)          # tilt of the C-X axis against the plane
    if hals and tilt > 15.0:
        print(f"  Note: the C-X axis of {hals[0]['label']} points "
              f"{tilt:.0f} degrees out of the best-fit plane;")
        print(f"    the sigma view follows the true bond axis, "
              f"pi/edge the plane.")
    # --- PyMOL scene ----------------------------------------------------
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
    cmd.disable("espramp")                 # keep the bar out of the image

    cmd.set("transparency", args.transparency)
    cmd.set("transparency_mode", 2)
    cmd.set("surface_quality", 1)
    cmd.set("two_sided_lighting", 1)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.15)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)              # no perspective -> comparable

    outdir = args.outdir or "."
    # Its own name suffix, otherwise a rainbow run overwrites the
    # red-white-blue image set of the same molecule.
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
            # Zoom on the molecule, NOT on "surf": the isosurface object
            # carries the extent of the whole grid box with it and would leave
            # the subject looking tiny.
            cmd.zoom("mol", BUFFER_ANGSTROM)

            suffix = f"_{bg}" if len(args.backgrounds) > 1 else ""
            png = os.path.join(outdir, f"{args.prefix}{cmap_tag}_{name}{suffix}.png")
            cmd.ray(args.width, args.height)
            cmd.png(png, dpi=args.dpi)
            written.append(png)
            print(f"    -> {png}")

    # --- colour bar as a separate image ---------------------------------
    bar = None
    try:
        bar = colorbar(os.path.join(outdir, f"{args.prefix}{cmap_tag}_colorbar.png"),
                       rng, dpi=args.dpi, rainbow=args.rainbow)
        written.append(bar)
        print(f"    -> {bar}")
    except ImportError:
        print("    (matplotlib is missing - the colour bar is skipped; "
              "'conda install matplotlib' to enable it)")

    # --- record ---------------------------------------------------------
    settings = os.path.join(outdir, f"{args.prefix}{cmap_tag}_settings.txt")
    with open(settings, "w", encoding="utf-8") as fh:
        fh.write("Render parameters (written by render_esp.py)\n")
        fh.write("=" * 55 + "\n")
        fh.write(f"Structure         : {args.struct}\n")
        fh.write(f"Density cube      : {args.density}\n")
        fh.write(f"ESP cube          : {args.esp}\n")
        fh.write(f"Grid              : {dens.shape[0]} x {dens.shape[1]} "
                 f"x {dens.shape[2]}\n")
        fh.write(f"Isovalue rho      : {args.iso} a.u.\n")
        fh.write(f"V_S,min           : {vmin:+.5f} a.u. "
                 f"({vmin*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))  on "
                 f"{loc.get('vmin_atom','?')}\n")
        fh.write(f"V_S,max           : {vmax:+.5f} a.u. "
                 f"({vmax*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))  on "
                 f"{loc.get('vmax_atom','?')}\n")
        # One line per halogen, sorted by sigma hole, descending.
        for e in hals:
            tag = f"({e['label']})"
            if "sigma_max" not in e:
                fh.write(f"sigma hole {tag:<7}: "
                         f"not evaluable (too few points)\n")
                continue
            fh.write(f"sigma hole {tag:<7}: "
                     f"{e['sigma_max']:+.5f} a.u. "
                     f"({e['sigma_max']*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))"
                     f"  [{e.get('sigma_method','point-based')}]\n")
            if "belt_min" in e:
                fh.write(f"belt       {tag:<7}: "
                         f"{e['belt_min']:+.5f} a.u. "
                         f"({e['belt_min']*HARTREE_TO_KCAL:+.2f} kcal/(mol*e))\n")
        if len(hals) > 1:
            fh.write(f"sigma view on     : {hals[0]['label']} "
                     f"(strongest sigma hole)\n")
        fh.write(f"Grid spacing      : {spacing:.4f} Bohr\n")
        fh.write(f"Colour scale      : -{rng:.4f} .. +{rng:.4f} a.u. ({how})\n")
        ramp_name = ("rainbow (red-yellow-green-cyan-blue)" if args.rainbow
                     else "red-white-blue")
        fh.write(f"Colour ramp       : {ramp_name}\n")
        fh.write(f"Transparency      : {args.transparency}\n")
        fh.write(f"Background        : {', '.join(args.backgrounds)}\n")
        fh.write(f"Image size        : {args.width} x {args.height} px, "
                 f"{args.dpi} dpi\n")
        fh.write(f"Projection        : orthoscopic\n")
        fh.write(f"Views             : {', '.join(views.keys())}\n")
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
        "sigma_method": loc.get("sigma_method", "point-based"),
        "sigma_angle": loc.get("sigma_angle"),
        "belt_min": loc.get("belt_min"),
        # All halogens, descending by sigma hole. The flat fields above refer
        # to halogens[0].
        "halogens": [{k: v for k, v in e.items() if k not in ("pos", "axis")}
                     for e in hals],
        "files": written,
        "settings_file": settings,
    }


def colorbar(path, rng, dpi=300, rainbow=False):
    """Horizontal colour bar as a separate PNG (needs matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    cmap = LinearSegmentedColormap.from_list(
        "esp", RAMP_HEX["rainbow" if rainbow else "redblue"])

    # A generous height plus bbox_inches="tight" when saving: otherwise the
    # axis labels at the bottom are cut off, which shows in the rendered
    # README.
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
        description="Render a standardised set of ESP images from cube files.")
    p.add_argument("--density", default=None, help="cube of the electron density")
    p.add_argument("--esp", default=None, help="cube of the ESP")
    p.add_argument("--struct", default=None,
                   help="structure file (.mol/.sdf/.xyz)")
    p.add_argument("--prefix", default=None, help="prefix of the image names")
    p.add_argument("--outdir", default="images", help="output folder")
    p.add_argument("--iso", type=float, default=0.001,
                   help="isovalue of the density surface in a.u. (default 0.001)")
    p.add_argument("--esp-range", default="auto",
                   help="'auto' or a fixed value in a.u., e.g. 0.03")
    p.add_argument("--transparency", type=float, default=0.15,
                   help="surface transparency 0..1 (default 0.15). "
                        "0 = opaque, clearest colours; 0.3+ makes the profile "
                        "and axial views unreadable, because you look through "
                        "the whole molecule.")
    p.add_argument("--backgrounds", nargs="+", default=["white"],
                   help="background colours, e.g. white black")
    p.add_argument("--width", type=int, default=2000)
    p.add_argument("--height", type=int, default=1600)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--rainbow", action="store_true",
                   help="rainbow colour ramp instead of red-white-blue. Red "
                        "stays negative, blue positive; yellow/green/cyan lie "
                        "in between. Writes an image set of its own, "
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
        raise SystemExit("Missing inputs: " + ", ".join(missing))

    if not args.prefix:
        base = os.path.splitext(os.path.basename(args.struct))[0]
        args.prefix = base.split("_")[0] or "molecule"

    print("=" * 70)
    print("render_esp.py - standardised ESP images")
    print("=" * 70)
    print(f"  structure : {args.struct}")
    print(f"  density   : {args.density}")
    print(f"  ESP       : {args.esp}")
    print(f"  prefix    : {args.prefix}")

    render_all(args)
    print("Done.")
    return 0


# Run as soon as the script is NOT imported as a module.
#
# Why not the usual  if __name__ == "__main__"  ?
# PyMOL executes .py files passed to it with exec() in a namespace of its own,
# in which __name__ is precisely not "__main__". With the standard check,
# nothing at all happens on  pymol -cq render_esp.py -- ... : the script is
# read, all functions are defined, and that is the end of it - without an error
# message. That silent non-execution is hard to diagnose, hence the inverted
# check here.
if __name__ != "render_esp":
    _argv = sys.argv[1:]
    if "--" in _argv:                  # called as: pymol -cq script.py -- ...
        _argv = _argv[_argv.index("--") + 1:]
    else:
        # On start-up PyMOL pushes arguments of its own into sys.argv. Throw
        # away everything before the script file, so that argparse does not
        # trip over it.
        for _i, _a in enumerate(_argv):
            if _a.endswith("render_esp.py"):
                _argv = _argv[_i + 1:]
                break
    main(_argv)
