#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xyzToCube.py
============

Converts Turbomole ``pointval`` grid files (td.xyz, tp.xyz, ...) into the
Gaussian cube format, so that they can be loaded into PyMOL, VMD, ChimeraX,
Avogadro or Multiwfn.

Background
----------
Turbomole writes volumetric data as a plain ASCII point cloud::

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

Every line carries the full coordinates -> at 251^3 points that is 1.25 GB. A
cube file stores the same information with an implicit grid (~200 MB).

Two pitfalls that this script catches:

1. **Axis order.** In the Turbomole file *x* varies fastest, in the cube format
   *z* varies fastest. Without reordering you get a transposed, mirrored
   molecule.
2. **Units.** The grid is in Bohr (atomic units), the structure file usually in
   Angstrom. By default the script converts the atoms
   (``--struct-unit angstrom``).

Usage
-----
Typical call for the bromobenzene example::

    python xyzToCube.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz --pymol

Result: ``td.cube``, ``tp.cube`` and ``esp.pml`` (a ready-to-run PyMOL script,
triggered by --pymol).

If PyMOL gets too slow with the full 251^3 grid, use every second point::

    python xyzToCube.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz --stride 2

``.xyz``, ``.mol`` and ``.sdf`` are accepted as structure files - the same
formats as in render_esp.py. Recommendation: give both scripts the *same*
file. Otherwise the atoms in the cube header come from one source and the
sticks in PyMOL from another, and a discrepancy between the two goes unnoticed
because no error is raised.

Only numpy is required (``pip install numpy``).
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
import time

import numpy as np

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

from constants import (BOHR_PER_ANGSTROM, ANGSTROM_PER_BOHR,  # noqa: F401
                       HARTREE_TO_KJ)

# Element symbols in order of atomic number: the position in the list IS Z-1,
# and SYMBOL_TO_Z is built from it below. The complete periodic table up to
# oganesson (Z = 118); the list used to end at radon, which would have created
# an avoidable failure case for the actinides.
ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYMBOL_TO_Z = {sym.upper(): i + 1 for i, sym in enumerate(ELEMENTS)}

CHUNK_BYTES = 1 << 24                      # 16 MB read buffer


# ----------------------------------------------------------------------------
# Reading the structure file
# ----------------------------------------------------------------------------

def _symbol_to_z(sym: str, path: str) -> int:
    """Element symbol (or atomic number) -> atomic number."""
    key = sym.strip().capitalize().upper()
    if key in SYMBOL_TO_Z:
        return SYMBOL_TO_Z[key]
    if re.fullmatch(r"\d+", sym.strip()):
        return int(sym)                                # atomic number, not a symbol
    raise ValueError(
        f"Unknown element '{sym}' in {path}. "
        f"Please extend the ELEMENTS list in the script."
    )


def _read_xyz(lines, path):
    """xyz format: lines of ``symbol x y z``.

    Accepts both the standard format (atom count + comment line + atoms) and a
    bare coordinate list, which Turbomole workflows often produce.
    """
    start = 0
    first = lines[0].strip() if lines else ""
    if re.fullmatch(r"\d+", first):
        start = 2                                     # skip count + comment

    atoms = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
        except ValueError:
            continue                                   # comment or junk line
        atoms.append((_symbol_to_z(parts[0], path), x, y, z))
    return atoms


def _read_molfile(lines, path):
    """MDL molfile / SD file (.mol, .sdf), V2000 and V3000.

    Layout of V2000::

        line 1    title
        line 2    program line
        line 3    comment
        line 4    counts line:  " 12 12  0 ... V2000"
        then      per atom:  x  y  z  symbol  ...     <- coordinates FIRST
        then      bond block

    For SD files only the first record is read (up to ``$$$$``).
    Molfile coordinates are in Angstrom by definition.
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
        raise ValueError(f"{path}: too short for a molfile.")

    counts = lines[3]
    try:
        natoms = int(counts[0:3])
    except ValueError:
        raise ValueError(
            f"{path}: counts line (line 4) not readable: {counts.strip()!r}")

    for line in lines[4:4 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: atom line incomplete: {line.strip()!r}")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        atoms.append((_symbol_to_z(parts[3], path), x, y, z))

    return atoms


def read_structure(path: str, unit: str = "angstrom"):
    """Reads a structure file and returns the atoms in **Bohr**.

    Supported:

    ==========  =====================================================
    ``.xyz``    ``symbol x y z``, with or without header lines
    ``.mol``    MDL molfile V2000/V3000 (coordinates *before* the symbol)
    ``.sdf``    SD file, first record
    ==========  =====================================================

    Returns a list of ``(Z, x, y, z)``.

    ``unit`` applies to xyz files only - molfile coordinates are in Angstrom by
    definition, and the setting is ignored there.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    if not lines:
        raise ValueError(f"{path} is empty.")

    ext = os.path.splitext(path)[1].lower()

    if ext in (".mol", ".sdf", ".sd"):
        # SD file: the first record only
        for i, ln in enumerate(lines):
            if ln.startswith("$$$$"):
                lines = lines[:i]
                break
        atoms = _read_molfile(lines, path)
        file_unit = "angstrom"
    else:
        atoms = _read_xyz(lines, path)
        file_unit = unit

    if not atoms:
        raise ValueError(f"No atoms found in {path}.")

    if file_unit == "angstrom":
        atoms = [(z, x * BOHR_PER_ANGSTROM,
                     y * BOHR_PER_ANGSTROM,
                     zz * BOHR_PER_ANGSTROM) for (z, x, y, zz) in atoms]
    elif file_unit != "bohr":
        raise ValueError("unit must be 'angstrom' or 'bohr'")

    return atoms


# ----------------------------------------------------------------------------
# Reading the Turbomole grid file
# ----------------------------------------------------------------------------

def parse_header(fh):
    """Reads the ``#`` header of a pointval file.

    Returns (info dict, first data line). The first data line has already been
    read from the stream and has to be processed by the caller.
    """
    info = {
        "origin": np.zeros(3),
        "vectors": np.eye(3),
        "grid": [None, None, None],       # each (start, delta, points)
        "title": "",
        "quantity": "",
    }
    first_data_line = None

    while True:
        line = fh.readline()
        if not line:
            raise ValueError("File ends inside the header - no data found.")
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
                raise ValueError(f"Grid line not readable: {line!r}")
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
        raise ValueError("Header incomplete: #grid1/#grid2/#grid3 are missing.")

    return info, first_data_line


def read_values(path, verbose=True):
    """Reads the 4th column of a pointval file as a float32 array.

    Reads in large blocks instead of line by line - for the 1.25 GB files that
    is about an order of magnitude faster.
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
            if cut == -1:                      # extremely long line - read on
                remainder = chunk
                continue
            remainder = chunk[cut + 1:]
            block = chunk[:cut]

            tokens = block.split()
            if not tokens:
                continue
            if len(tokens) % 4 != 0:
                raise ValueError(
                    f"{path}: expected 4 columns per line, "
                    f"found {len(tokens)} values in one block."
                )
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            n = arr.size
            if filled + n > total:
                raise ValueError(
                    f"{path}: more data points than the header states "
                    f"({filled + n} > {total})."
                )
            values[filled:filled + n] = arr
            filled += n

            if verbose:
                pct = 100.0 * filled / total
                sys.stdout.write(f"\r    reading {os.path.basename(path)}: "
                                 f"{pct:5.1f} %")
                sys.stdout.flush()

        # the last, incompletely buffered line
        tokens = remainder.split()
        if tokens:
            if len(tokens) % 4 != 0:
                raise ValueError(f"{path}: last line incomplete.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            values[filled:filled + arr.size] = arr
            filled += arr.size

    if verbose:
        sys.stdout.write(f"\r    reading {os.path.basename(path)}: 100.0 %  "
                         f"({filled:,} points in {time.time() - t0:.1f} s)\n")

    if filled != total:
        raise ValueError(
            f"{path}: read {filled} values, header states {total}."
        )

    # Turbomole: x varies fastest -> the memory layout is [i3, i2, i1]
    data = values.reshape(n3, n2, n1)
    # Cube: z varies fastest -> we want [i1, i2, i3]
    data = np.ascontiguousarray(np.transpose(data, (2, 1, 0)))

    return info, data


# ----------------------------------------------------------------------------
# Writing the cube
# ----------------------------------------------------------------------------

def write_cube(path, info, data, atoms, stride=1, comment=""):
    """Writes a Gaussian cube. All lengths in Bohr.

    It is written to <name>.part first and renamed at the end. A full cube is
    200 MB and takes minutes; if the run is aborted in that time - Ctrl-C, a
    closed window, a full disk - half a file with a valid header would
    otherwise be left behind. The next run takes that for finished, skips the
    conversion and renders an isosurface with its rear half missing. The
    rename is the moment the file comes into being - before that it does not
    exist under its name.
    """
    if stride > 1:
        data = data[::stride, ::stride, ::stride]

    n = data.shape
    starts = np.array([info["grid"][i][0] for i in range(3)])
    deltas = np.array([info["grid"][i][1] for i in range(3)])
    vecs = info["vectors"]

    # Ursprung des ersten Voxels im kartesischen Raum
    origin = info["origin"] + sum(starts[i] * vecs[i] for i in range(3))
    # voxel vectors (scaled by the stride)
    voxel = np.array([deltas[i] * stride * vecs[i] for i in range(3)])

    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(f"{comment or 'cube written by xyzToCube.py'}\n")
        fh.write(f"{info.get('quantity', '') or 'volumetric data'} | "
                 f"{info.get('title', '')} | units: Bohr\n")
        fh.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} "
                 f"{origin[2]:12.6f}\n")
        for i in range(3):
            fh.write(f"{n[i]:5d} {voxel[i][0]:12.6f} {voxel[i][1]:12.6f} "
                     f"{voxel[i][2]:12.6f}\n")
        for (znum, x, y, z) in atoms:
            fh.write(f"{znum:5d} {float(znum):12.6f} {x:12.6f} {y:12.6f} "
                     f"{z:12.6f}\n")

        # Values: z fastest, 6 per line.
        # A precompiled format pattern per z row is considerably faster than a
        # loop over all 15.8 million individual values.
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

    # Only here does <name>.cube exist. os.replace also replaces an existing
    # file and is atomic within one file system - there is no moment in which
    # the old cube is gone and the new one is not yet there.
    os.replace(tmp, path)
    return n, origin, voxel


# ----------------------------------------------------------------------------
# Shell evaluation (rho = iso) and colour scale
#
# The single source for both scripts: render_esp.py imports shell_mask(),
# esp_statistics() and nice_range() from here instead of defining them a second
# time. The mask logic used to sit in the repository three times (here in
# auto_esp_range, in render_esp.esp_statistics and in render_esp.shell_points).
# A changed shell thickness would have had to be followed through in all three
# places - otherwise the numbers in the statistics and the colour scale in the
# picture would have drifted apart silently.
# ----------------------------------------------------------------------------

# Shell thickness relative to the isovalue. 0.12 keeps the shell thin enough
# that the values really come from the isosurface; if that leaves it too
# sparsely populated (coarse grid, small molecule), SHELL_TOL_FALLBACK widens it
# once.
SHELL_TOL_FACTOR = 0.12
SHELL_TOL_FALLBACK = 0.30
SHELL_MIN_POINTS = 50


def shell_mask(density, iso=0.001, tol_factor=SHELL_TOL_FACTOR):
    """Boolean mask of the grid points on the rho = iso shell.

    Only the mask is returned, not the values - that way callers can apply the
    same selection to ESP values, grid indices or coordinates without copying
    the grid, which at 251^3 is several hundred MB.
    """
    mask = np.abs(density - iso) < iso * tol_factor
    if mask.sum() < SHELL_MIN_POINTS:         # shell too thin -> widen it
        mask = np.abs(density - iso) < iso * SHELL_TOL_FALLBACK
    return mask


def esp_statistics(density, esp, iso=0.001, tol_factor=SHELL_TOL_FACTOR):
    """ESP statistics on the rho = iso shell.

    Returns (V_min, V_max, number_of_points) in atomic units.
    ``tol_factor`` sets the shell thickness relative to the isovalue.
    """
    mask = shell_mask(density, iso, tol_factor)
    if not mask.any():
        return None, None, 0
    shell = esp[mask]
    return float(shell.min()), float(shell.max()), int(mask.sum())


def nice_range(vmin, vmax, step=0.005):
    """A symmetric range, rounded up to a multiple of ``step``."""
    amp = max(abs(vmin), abs(vmax))
    return math.ceil(amp / step) * step


def auto_esp_range(density, esp, iso=0.001, tol_factor=SHELL_TOL_FACTOR,
                   step=0.005):
    """A symmetric colour scale from the ESP values on the rho = iso shell.

    Do NOT take the range from the whole grid (there the nuclear singularities
    dominate with several hundred a.u.), take it only from the points on the
    isosurface, and round the result up symmetrically to a round multiple of
    ``step``. Exactly the combination render_esp.py uses for
    ``--esp-range auto`` - from the same two functions over there.

    Returns the half width in a.u., or None if no shell is found (then the
    density file is missing, or the isovalue does not match the data).
    """
    if density is None or esp is None or density.shape != esp.shape:
        return None
    vmin, vmax, npts = esp_statistics(density, esp, iso, tol_factor)
    if npts == 0:
        return None
    return nice_range(vmin, vmax, step)


# ----------------------------------------------------------------------------
# Does the structure belong to the grid?
# ----------------------------------------------------------------------------

class StructureGridMismatch(Exception):
    """The structure file is not the geometry the grids were computed from."""


# Median electron density at a heavy nucleus, below which the two do not
# belong together. Measured across the molecules of this project:
#
#     grid spacing 0.12 Bohr   median  64
#     grid spacing 0.25 Bohr   median  27
#     grid spacing 0.60 Bohr   median   5      the decimated reference set
#     ------------------------------------
#     structure from a different run       0.03
#
# A coarse grid pushes the value down, because the nucleus then sits further
# from the nearest grid point and the density peak falls off steeply. 1.0 is a
# factor of five below the coarsest healthy case and thirty above the broken
# one; there is nothing in between that could be mistaken for the other.
NUCLEI_RHO_MIN = 1.0


def nuclei_density(density, atoms, origin, voxel):
    """Median electron density at the heavy nuclei, or None if not applicable.

    A cube file marries two things from different sources: the atom list comes
    from the structure file, the grid from the quantum chemistry. Nothing in
    the format forces them to agree, and when they do not, the molecule floats
    beside its own isosurface - and every number of the run refers to atoms
    that are not where the density says they are.

    The check is one line of physics: where a nucleus sits, the electron
    density is enormous (about 60 for carbon on a fine grid), and where there
    is none it is near zero. So look up the density at every atom position; if
    the two belong together, the nuclei sit on peaks.

    Two details matter:

    * The MEDIAN, not the minimum. Iodine is described by an effective core
      potential, so its core electrons are not in the density at all and its
      nucleus carries about 0.26 - less than a hydrogen, on a perfectly
      healthy molecule. A minimum would reject every iodine compound. The
      median tolerates a minority of such atoms.
    * The largest value in the 3x3x3 neighbourhood, not the value at the
      nearest point. A nucleus rarely sits exactly on a grid point, and this
      removes the luck of where it falls. Trilinear interpolation would be the
      wrong tool here: the density has a cusp at the nucleus, and interpolation
      cuts the peak off - measured at 46 % too low on the coarse grid.
    """
    if density is None or not atoms:
        return None
    vox = np.asarray(voxel, dtype=float)
    if vox.ndim == 2:
        diag = np.diag(vox)
        if not np.allclose(vox, np.diag(diag)):
            return None                    # not axis-aligned - no cheap lookup
        vox = diag
    shape = np.array(density.shape)
    vals = []
    for atom in atoms:
        z = atom[0]
        if z <= 1:
            continue                       # hydrogen carries too little
        i = np.rint((np.array(atom[1:4]) - origin) / vox).astype(int)
        if (i < 0).any() or (i >= shape).any():
            vals.append(0.0)               # outside the grid counts as a miss
            continue
        lo = np.maximum(i - 1, 0)
        hi = np.minimum(i + 2, shape)
        vals.append(float(density[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]].max()))
    return float(np.median(vals)) if vals else None


def check_alignment(density, atoms, origin, voxel, label=""):
    """Raise StructureGridMismatch if structure and grid do not belong together."""
    med = nuclei_density(density, atoms, origin, voxel)
    if med is None or med >= NUCLEI_RHO_MIN:
        return med
    where = f"{label}: " if label else ""
    raise StructureGridMismatch(
        f"{where}structure and grid do not match. The median electron density "
        f"at the heavy nuclei is {med:.3f}, expected > {NUCLEI_RHO_MIN:g} - "
        f"the atoms sit in near-empty space.\n"
        f"    The structure file is probably not the geometry the grids were "
        f"computed from (a different conformer, a different orientation, or a "
        f"left-over file from an earlier run).")

# ----------------------------------------------------------------------------
# Generating the PyMOL script
# ----------------------------------------------------------------------------

PML_TEMPLATE = """# --------------------------------------------------------------
# esp.pml - ESP on the electron density isosurface
# Start:  pymol esp.pml       (or inside PyMOL:  @esp.pml)
# --------------------------------------------------------------

reinitialize

# 1) load the structure and the volumetric data
load {struct}, mol
{load_density}
load {esp_cube}, esp

# 2) the molecule as sticks
hide everything
show sticks, mol
set stick_radius, 0.3
color grey70, mol and elem C
util.cnc mol

# 3) isosurface of the electron density at rho = {iso} a.u.
#    (the Politzer/Murray convention for the "molecular surface")
{isosurface}

# 4) colour ramp for the ESP; values in Hartree/e (a.u.)
#    {vmin} .. {vmax} a.u.  equals {kvmin:.0f} .. {kvmax:.0f} kJ/(mol*e)
ramp_new espramp, esp, [{ramp_levels}], [{ramp_colors}]

# 5) map the ESP onto the surface
set surface_color, espramp, {surface_target}
set surface_quality, 1

#    transparency: 0 = opaque (strongest colours, sticks invisible),
#    0.15 = the default (the skeleton shows through),
#    from about 0.3 on it becomes unreadable, because you look through the
#    whole molecule.
set transparency, {transparency}
set transparency_mode, 2
set two_sided_lighting, on

# 6) appearance / rendering
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 0
set specular, 0.2
set ambient, 0.15
orient mol
zoom mol, 2.0

# 7) high-resolution image
# ray 2400, 1800
# png esp.png, dpi=300
"""


# Colour ramps - identical to render_esp.RAMP_PYMOL. Red stays negative in
# both and blue positive; the rainbow only pushes yellow/green/cyan in between,
# so that the two image sets can be laid side by side.
PML_RAMPS = {
    "redblue": ["red", "white", "blue"],
    "rainbow": ["red", "yellow", "green", "cyan", "blue"],
}


def write_pymol_script(path, struct, density_cube, esp_cube, vmin, vmax,
                       iso=0.001, transparency=0.15, rainbow=False):
    if density_cube:
        load_density = f"load {density_cube}, dens"
        isosurface = f"isosurface surf, dens, {iso}"
        surface_target = "surf"
    else:
        load_density = "# no density file available"
        isosurface = ("# Substitute: van der Waals surface from the structure.\n"
                      "# The Politzer/Murray convention (rho = 0.001 a.u.)\n"
                      "# would need the density file td.cube.\n"
                      "show surface, mol\n"
                      "set surface_solvent, 0")
        surface_target = "mol"

    cols = PML_RAMPS["rainbow" if rainbow else "redblue"]

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
        kvmin=vmin * HARTREE_TO_KJ, kvmax=vmax * HARTREE_TO_KJ,
        ramp_levels=", ".join(
            f"{vmin + (vmax - vmin) * i / (len(cols) - 1):g}"
            for i in range(len(cols))),
        ramp_colors=", ".join(cols),
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


# ----------------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Convert Turbomole pointval grid files (td.xyz, tp.xyz) "
                    "to Gaussian cube.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
               "  python xyzToCube.py --struct brombenzol_aro_opti.mol "
               "td.xyz tp.xyz --pymol\n",
    )
    p.add_argument("grids", nargs="+",
                   help="Turbomole grid files (e.g. td.xyz tp.xyz)")
    p.add_argument("--struct", "-s", required=True,
                   help="structure file: .xyz, .mol or .sdf")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom",
                   help="unit of the structure file (default: angstrom); applies to .xyz only - .mol/.sdf are always Angstrom")
    p.add_argument("--outdir", "-o", default=None,
                   help="output directory (default: next to the input)")
    p.add_argument("--stride", type=int, default=1,
                   help="write every n-th grid point only "
                        "(2 => 8x smaller file, default: 1)")
    p.add_argument("--quiet", "-q", action="store_true")

    # These options change NOTHING about the cube files - they only describe
    # the scene that --pymol writes on the side. Hence a group of their own and
    # the --pml- prefix on everything that could otherwise be confused with an
    # option of the same name in render_esp.py.
    g = p.add_argument_group("PyMOL scene (only together with --pymol)")
    g.add_argument("--pymol", action="store_true",
                   help="also write a ready-to-run esp.pml")
    g.add_argument("--esp-range", default="auto",
                   help="half width of the ESP colour scale in a.u., or 'auto' "
                        "(default): determined from the ESP values on the "
                        "isosurface, as in render_esp.py")
    g.add_argument("--pml-iso", type=float, default=0.001,
                   help="isovalue of the density surface INSIDE esp.pml "
                        "(default: 0.001). Changes the scene only, not the "
                        "cube data - unlike --iso in render_esp.py, which "
                        "shifts the measured statistics.")
    g.add_argument("--transparency", type=float, default=0.15,
                   help="surface transparency in esp.pml, 0..1 "
                        "(default: 0.15; 0 = opaque)")
    g.add_argument("--rainbow", action="store_true",
                   help="rainbow ramp in esp.pml instead of red-white-blue")
    args = p.parse_args(argv)

    verbose = not args.quiet

    if verbose:
        print("=" * 70)
        print("xyzToCube.py - Turbomole pointval  ->  Gaussian Cube")
        print("=" * 70)

    atoms = read_structure(args.struct, unit=args.struct_unit)
    if verbose:
        print(f"[1] structure: {args.struct} -> {len(atoms)} atoms "
              f"(read as {args.struct_unit}, stored as Bohr)")

    # For --esp-range auto the density and the ESP are needed once more after
    # writing. Keep them in memory only then: at 251^3 that is 63 MB per grid,
    # which would otherwise sit around for nothing.
    need_auto = args.pymol and str(args.esp_range).lower() == "auto"
    cache = {}

    written = {}
    for gpath in args.grids:
        if verbose:
            print(f"[2] grid file: {gpath}")
        info, data = read_values(gpath, verbose=verbose)

        n1, n2, n3 = (info["grid"][i][2] for i in range(3))
        if verbose:
            print(f"    grid {n1} x {n2} x {n3}, "
                  f"delta = {info['grid'][0][1]} Bohr, "
                  f"quantity = '{info['quantity'] or 'unknown'}'")
            print(f"    value range: {data.min():+.6g} .. {data.max():+.6g}")

        base = os.path.splitext(os.path.basename(gpath))[0]
        outdir = args.outdir or os.path.dirname(os.path.abspath(gpath))
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, base + ".cube")

        shape, origin, voxel = write_cube(
            outpath, info, data, atoms, stride=args.stride,
            comment=f"{info['quantity'] or base} - converted from {os.path.basename(gpath)}",
        )
        if verbose:
            mb = os.path.getsize(outpath) / 1024 ** 2
            print(f"    -> {outpath}  ({shape[0]}x{shape[1]}x{shape[2]}, "
                  f"{mb:.1f} MB)")

        q = (info["quantity"] or "").lower()
        if "potential" in q:
            written["esp"] = outpath
            if need_auto:
                cache["esp"] = data
        elif "density" in q:
            written["density"] = outpath
            if need_auto:
                cache["density"] = data
        else:
            written.setdefault("other", []).append(outpath)

    if args.pymol:
        outdir = args.outdir or os.path.dirname(os.path.abspath(args.grids[0]))
        pml = os.path.join(outdir, "esp.pml")
        esp_cube = written.get("esp")
        if esp_cube is None:
            print("    ! no ESP file recognised - esp.pml is skipped.",
                  file=sys.stderr)
        else:
            rng = None
            if str(args.esp_range).lower() == "auto":
                rng = auto_esp_range(cache.get("density"), cache.get("esp"),
                                     iso=args.pml_iso)
                if rng is None:
                    rng = 0.03
                    print("    ! colour scale not automatically determinable "
                          "(no density file?) - using +/- 0.03 a.u.",
                          file=sys.stderr)
                elif verbose:
                    print(f"    colour scale, automatic: +/- {rng:.4f} a.u.")
            else:
                rng = float(args.esp_range)

            write_pymol_script(
                pml,
                struct=os.path.relpath(os.path.abspath(args.struct), outdir),
                density_cube=(os.path.basename(written["density"])
                              if "density" in written else None),
                esp_cube=os.path.basename(esp_cube),
                vmin=-rng, vmax=rng, iso=args.pml_iso,
                transparency=args.transparency, rainbow=args.rainbow,
            )
            if verbose:
                print(f"[3] PyMOL script: {pml}")
                print(f"    start with:  pymol {os.path.basename(pml)}")

    if verbose:
        print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
