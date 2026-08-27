#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CreateTpTdFromSmiles.py
=======================

Produces test datasets for the ESP workflow from a SMILES string.

The output per molecule is a folder holding exactly the files that come out of
a Turbomole run::

    <name>/
        td.xyz        electron density, Turbomole pointval format
        tp.xyz        electrostatic potential, the same format
        <name>.mol    structure (MDL molfile, Angstrom)

The folder can be dropped straight into ``sandbox/``; ``run_all.py`` finds it
there and the pipeline runs through unchanged from ``xyzToCube.py`` on. That is
deliberate: if cube files were written here directly, the converter with its
unit conversion and index reordering would stay untested - precisely the part
most likely to break.

Flow::

    SMILES --RDKit--> 3D geometry --PySCF--> density matrix --> rho(r), V(r)

Details, choice of parameters and limitations are in ``README.txt`` in the same
folder. In short: these are **test fixtures, not reference data**.

Call
----
    python CreateTpTdFromSmiles.py --preset
    python CreateTpTdFromSmiles.py --smiles "CC(=O)c1ccc(Br)cc1" --name bromacetophenon
    python CreateTpTdFromSmiles.py --preset --spacing 0.35 --margin 3.0

Dependencies: rdkit, pyscf, numpy (see environment-testdata.yml).
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

# Use the same conversion factor as the pipeline in scripts/. Were there two
# values, the atoms of the generated .mol would sit slightly differently from
# the generated grid - invisible in the picture, wrong in the numbers.
# scripts/constants deliberately has no dependencies and can therefore be
# imported from the separate esp-testdata environment as well.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "scripts"))
from constants import BOHR_PER_ANGSTROM  # noqa: E402

# Built-in test cases. Chosen deliberately to cover the gaps that three nearly
# identical halobenzenes leave open.
PRESETS = [
    # Halogen AND carbonyl: V_S,min has to move to the carbonyl oxygen while
    # the belt value stays on the bromine. Only here do the two numbers diverge
    # at all.
    ("4-bromacetophenon", "CC(=O)c1ccc(Br)cc1"),
    # No halogen: does the sigma-hole analysis run into nothing cleanly, and
    # does the orientation fall back sensibly onto the principal axes without a
    # C-X axis?
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
]


# ----------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------

def build_geometry(smiles, name, outdir, seed=0xF00D):
    """SMILES -> 3D geometry -> molfile. Returns (atom list, molfile path).

    The coordinates come from ETKDG plus MMFF94 optimisation, that is, from a
    force field. That is enough for a functional test of the pipeline; not for
    numbers meant to stand beside quantum-chemically optimised geometries.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"SMILES not readable: {smiles!r}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed                 # reproducible conformation
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise SystemExit(f"{name}: 3D embedding failed.")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        print("    ! MMFF does not know all atom types, using UFF")

    os.makedirs(outdir, exist_ok=True)
    molfile = os.path.join(outdir, f"{name}.mol")
    Chem.MolToMolFile(mol, molfile)

    conf = mol.GetConformer()
    atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        atoms.append((atom.GetSymbol(), p.x, p.y, p.z))   # Angstrom
    return atoms, molfile


# ----------------------------------------------------------------------------
# Quantum chemistry
# ----------------------------------------------------------------------------

def run_scf(atoms, basis="def2-svp", method="hf", verbose=0):
    """SCF calculation. Returns (Mole object, density matrix).

    What we are after is not the energy but the density matrix - rho(r) and
    V(r) are evaluated on the grid from it afterwards.
    """
    from pyscf import gto, scf, dft

    atom_spec = [(s, (x, y, z)) for (s, x, y, z) in atoms]
    mol = gto.M(atom=atom_spec, basis=basis, unit="Angstrom", verbose=verbose)

    if method.lower() in ("b3lyp", "dft"):
        mf = dft.RKS(mol)
        mf.xc = "b3lyp"
    else:
        mf = scf.RHF(mol)
    mf.conv_tol = 1e-9
    energy = mf.kernel()
    if not mf.converged:
        print("    ! SCF did not converge - treat the result with caution")
    return mol, mf.make_rdm1(), energy


# ----------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------

def make_grid(atoms, spacing_bohr=0.25, margin_ang=3.5):
    """A regular grid around the molecule, everything in Bohr.

    The margin has to enclose the rho = 0.001 isosurface safely; that sits at
    roughly the van der Waals radius, so a good 2 Angstrom beyond the outermost
    nuclei. A margin of 3.5 Angstrom leaves enough room without inflating the
    grid needlessly - every additional Angstrom costs quadratically in
    computing time for the potential evaluation.
    """
    coords = np.array([[x, y, z] for (_, x, y, z) in atoms]) * BOHR_PER_ANGSTROM
    margin = margin_ang * BOHR_PER_ANGSTROM

    lo = coords.min(axis=0) - margin
    hi = coords.max(axis=0) + margin
    npts = np.ceil((hi - lo) / spacing_bohr).astype(int) + 1
    return lo, spacing_bohr, npts


def evaluate(mol, dm, lo, delta, npts, chunk=None, mem_mb=400, label=""):
    """Evaluates rho(r) and V(r) on the grid.

    V(r) is computed block by block via ``int1e_grids`` - the PySCF routine
    meant for exactly this, potential evaluation at many points. Block by
    block, because otherwise the intermediate matrix would be the memory hog:
    one full basis-function matrix per grid point.
    """
    from pyscf.dft import numint

    # Derive the block size from the memory need, do not guess it.
    # int1e_grids returns an array of shape (points, nao, nao) - at 193 basis
    # functions that is already 0.28 MB PER grid point. With a fixed block of
    # 20000 points this function tried to allocate 6 GB and was killed by the
    # operating system. Hence: choose the block so that the intermediate matrix
    # stays below mem_mb.
    nao = mol.nao
    if chunk is None:
        chunk = max(64, int(mem_mb * 1024**2 / (nao * nao * 8)))

    ax = [lo[i] + np.arange(npts[i]) * delta for i in range(3)]
    # Order as in Turbomole: x varies fastest, z slowest
    Z, Y, X = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")
    grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)
    total = grid.shape[0]

    rho = np.empty(total)
    esp = np.empty(total)

    charges = mol.atom_charges().astype(float)
    nuclei = mol.atom_coords()                       # Bohr

    if label:
        print(f"    {nao} basis functions -> block of {chunk:,} points "
              f"({chunk * nao * nao * 8 / 1024**2:.0f} MB intermediate matrix)")

    t0 = time.time()
    for start in range(0, total, chunk):
        pts = grid[start:start + chunk]

        # electron density
        ao = numint.eval_ao(mol, pts)
        rho[start:start + chunk] = numint.eval_rho(mol, ao, dm)

        # electrostatic potential: nuclei minus electrons
        d = np.linalg.norm(pts[:, None, :] - nuclei[None, :, :], axis=2)
        d[d < 1e-12] = 1e-12
        v_nuc = (charges[None, :] / d).sum(axis=1)

        # int1e_grids returns <i|1/|r-C||j> with a POSITIVE sign, so the
        # electronic contribution to the potential has to be subtracted.
        # Check: far outside a neutral molecule V has to go to 0. For HCl at a
        # distance of 40 Bohr: v_nuc = +0.4516, v_ele = +0.4512, difference
        # +0.0004 - the sum would be +0.9028 and thus obviously wrong.
        ints = mol.intor("int1e_grids", grids=pts)
        v_ele = np.einsum("pij,ij->p", ints, dm)

        esp[start:start + chunk] = v_nuc - v_ele

        done = min(start + chunk, total)
        if label:
            pct = 100.0 * done / total
            eta = (time.time() - t0) / done * (total - done)
            sys.stdout.write(f"\r    {label}: {pct:5.1f} %  "
                             f"(~{eta / 60:.1f} min left)   ")
            sys.stdout.flush()
    if label:
        sys.stdout.write(f"\r    {label}: 100.0 %  "
                         f"({total:,} points in {(time.time()-t0)/60:.1f} min)\n")

    return rho, esp


# ----------------------------------------------------------------------------
# Writing the Turbomole pointval format
# ----------------------------------------------------------------------------

def write_pointval(path, values, lo, delta, npts, quantity, title):
    """Writes a grid in the Turbomole ``pointval`` format.

    Deliberately the same format as the real data, including the full
    coordinates per line and the x-fastest order. Only that way does the test
    case pass through xyzToCube.py as well.
    """
    ax = [lo[i] + np.arange(npts[i]) * delta for i in range(3)]
    Z, Y, X = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")

    header = [
        f"#origin        {0.0:14.6f}{0.0:14.6f}{0.0:14.6f}",
        f"#vector1       {1.0:14.6f}{0.0:14.6f}{0.0:14.6f}",
        f"#vector2       {0.0:14.6f}{1.0:14.6f}{0.0:14.6f}",
        f"#vector3       {0.0:14.6f}{0.0:14.6f}{1.0:14.6f}",
    ]
    for i in range(3):
        header.append(f"#grid{i+1}  start  {lo[i]:.6f}  delta    {delta:.6f}"
                      f"  points    {npts[i]}")
    header += [
        f"#title for this grid {title}",
        f"#{quantity}",
        "#plotdata",
        "# cartesian coordinates x,y,z and f(x,y,z)",
    ]

    cols = np.stack([X.ravel(), Y.ravel(), Z.ravel(), values], axis=-1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(header) + "\n")
        fmt = "%20.8f%20.8f%20.8f%20.8f\n"
        buf = []
        for row in cols:
            buf.append(fmt % tuple(row))
            if len(buf) >= 8192:
                fh.write("".join(buf))
                buf.clear()
        if buf:
            fh.write("".join(buf))


# ----------------------------------------------------------------------------
# One molecule, end to end
# ----------------------------------------------------------------------------

def make_case(name, smiles, root, spacing, margin, basis, method):
    print(f"\n[{name}]  {smiles}")
    outdir = os.path.join(root, name)

    atoms, molfile = build_geometry(smiles, name, outdir)
    print(f"    geometry: {len(atoms)} atoms -> {os.path.basename(molfile)}")

    t0 = time.time()
    mol, dm, energy = run_scf(atoms, basis=basis, method=method)
    print(f"    SCF ({method}/{basis}): E = {energy:.6f} Hartree, "
          f"{mol.nao} basis functions, {time.time()-t0:.1f} s")

    lo, delta, npts = make_grid(atoms, spacing_bohr=spacing, margin_ang=margin)
    total = int(np.prod(npts))
    print(f"    grid: {npts[0]} x {npts[1]} x {npts[2]} = {total:,} points, "
          f"delta = {delta} Bohr")

    rho, esp = evaluate(mol, dm, lo, delta, npts, label="rho and V")

    td = os.path.join(outdir, "td.xyz")
    tp = os.path.join(outdir, "tp.xyz")
    write_pointval(td, rho, lo, delta, npts, "density", "101")
    write_pointval(tp, esp, lo, delta, npts, "electrostatic potential", "111")

    for p in (td, tp):
        print(f"    -> {p}  ({os.path.getsize(p)/1024**2:.0f} MB)")
    print(f"    rho: {rho.min():.4g} .. {rho.max():.4g}   "
          f"V: {esp.min():+.4g} .. {esp.max():+.4g}  (a.u.)")
    return outdir


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Produce test datasets (td.xyz/tp.xyz) from SMILES.")
    p.add_argument("--smiles", help="SMILES of a single molecule")
    p.add_argument("--name", help="folder and file name for it")
    p.add_argument("--preset", action="store_true",
                   help="produce the two built-in test cases")
    p.add_argument("--outdir", default="../sandbox",
                   help="root directory for the molecule folders")
    p.add_argument("--spacing", type=float, default=0.25,
                   help="grid spacing in Bohr (default 0.25)")
    p.add_argument("--margin", type=float, default=3.5,
                   help="margin around the molecule in Angstrom (default 3.5)")
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--method", default="hf", choices=["hf", "b3lyp"])
    args = p.parse_args(argv)

    print("=" * 70)
    print("CreateTpTdFromSmiles.py - test data for the ESP workflow")
    print("=" * 70)
    print("Note: force-field geometry, its own level of theory.")
    print("These are test fixtures - not reference data for comparison tables.")

    cases = list(PRESETS) if args.preset else []
    if args.smiles:
        if not args.name:
            raise SystemExit("--smiles also needs --name")
        cases.append((args.name, args.smiles))
    if not cases:
        raise SystemExit("Nothing to do: pass --preset or --smiles/--name.")

    made = []
    for name, smiles in cases:
        made.append(make_case(name, smiles, args.outdir, args.spacing,
                              args.margin, args.basis, args.method))

    print("\nDone. Folders produced:")
    for m in made:
        print(f"  {m}")
    print("\nContinue with:")
    print("  cd ../scripts && python run_all.py --root ../sandbox --stride 1")
    return 0


if __name__ != "CreateTpTdFromSmiles":
    _argv = sys.argv[1:]
    for _i, _a in enumerate(_argv):
        if _a.endswith("CreateTpTdFromSmiles.py"):
            _argv = _argv[_i + 1:]
            break
    main(_argv)
