#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CreateTpTdFromSmiles.py
=======================

Erzeugt Testdatensaetze fuer den ESP-Workflow aus einem SMILES-String.

Ausgabe pro Molekuel ist ein Ordner mit genau den Dateien, die auch aus einem
Turbomole-Lauf kommen::

    <name>/
        td.xyz        Elektronendichte, Turbomole-pointval-Format
        tp.xyz        elektrostatisches Potential, dasselbe Format
        <name>.mol    Struktur (MDL-Molfile, Angstrom)

Der Ordner laesst sich direkt nach ``sandbox/`` legen; ``run_all.py`` findet ihn
dort und die Pipeline laeuft ab ``xyzToCube.py`` unveraendert durch. Das ist
Absicht: wuerden hier gleich Cube-Dateien geschrieben, bliebe der Konverter
samt Einheitenumrechnung und Index-Umsortierung ungetestet - also genau der
Teil, der am ehesten bricht.

Ablauf::

    SMILES --RDKit--> 3D-Geometrie --PySCF--> Dichtematrix --> rho(r), V(r)

Details, Parameterwahl und Einschraenkungen stehen in ``README.txt`` im selben
Ordner. Kurzfassung: das sind **Testfixtures, keine Referenzdaten**.

Aufruf
------
    python CreateTpTdFromSmiles.py --preset
    python CreateTpTdFromSmiles.py --smiles "CC(=O)c1ccc(Br)cc1" --name bromacetophenon
    python CreateTpTdFromSmiles.py --preset --spacing 0.35 --margin 3.0

Abhaengigkeiten: rdkit, pyscf, numpy (siehe environment-testdata.yml).
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

BOHR_PER_ANGSTROM = 1.8897259886

# Eingebaute Testfaelle. Bewusst so gewaehlt, dass sie die Luecken abdecken,
# die drei fast identische Halogenbenzole offenlassen.
PRESETS = [
    # Halogen UND Carbonyl: V_S,min muss auf den Carbonyl-Sauerstoff wandern,
    # waehrend der Guertelwert am Brom bleibt. Erst hier laufen die beiden
    # Kennwerte ueberhaupt auseinander.
    ("4-bromacetophenon", "CC(=O)c1ccc(Br)cc1"),
    # Kein Halogen: laeuft die sigma-Loch-Analyse sauber ins Leere, und faellt
    # die Orientierung ohne C-X-Achse vernuenftig auf die Hauptachsen zurueck?
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
]


# ----------------------------------------------------------------------------
# Geometrie
# ----------------------------------------------------------------------------

def build_geometry(smiles, name, outdir, seed=0xF00D):
    """SMILES -> 3D-Geometrie -> Molfile. Rueckgabe: (atomliste, molfile-pfad).

    Die Koordinaten stammen aus ETKDG plus MMFF94-Optimierung, also aus einem
    Kraftfeld. Fuer einen Funktionstest der Pipeline reicht das; fuer Zahlen,
    die neben quantenchemisch optimierten Geometrien stehen sollen, nicht.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"SMILES nicht lesbar: {smiles!r}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed                 # reproduzierbare Konformation
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise SystemExit(f"{name}: 3D-Einbettung fehlgeschlagen.")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        print("    ! MMFF kennt nicht alle Atomtypen, UFF benutzt")

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
# Quantenchemie
# ----------------------------------------------------------------------------

def run_scf(atoms, basis="def2-svp", method="hf", verbose=0):
    """SCF-Rechnung. Rueckgabe: (Mole-Objekt, Dichtematrix).

    Uns interessiert nicht die Energie, sondern die Dichtematrix - daraus
    werden anschliessend rho(r) und V(r) auf dem Gitter ausgewertet.
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
        print("    ! SCF nicht konvergiert - Ergebnis mit Vorsicht behandeln")
    return mol, mf.make_rdm1(), energy


# ----------------------------------------------------------------------------
# Gitter
# ----------------------------------------------------------------------------

def make_grid(atoms, spacing_bohr=0.25, margin_ang=3.5):
    """Regelmaessiges Gitter um das Molekuel, alles in Bohr.

    Der Rand muss die rho=0.001-Isoflaeche sicher einschliessen; die liegt
    ungefaehr auf dem van-der-Waals-Radius, also gut 2 Angstrom jenseits der
    aeussersten Kerne. 3.5 Angstrom Rand lassen genug Luft, ohne das Gitter
    unnoetig aufzublaehen - jeder zusaetzliche Angstrom kostet quadratisch
    Rechenzeit bei der Potentialauswertung.
    """
    coords = np.array([[x, y, z] for (_, x, y, z) in atoms]) * BOHR_PER_ANGSTROM
    margin = margin_ang * BOHR_PER_ANGSTROM

    lo = coords.min(axis=0) - margin
    hi = coords.max(axis=0) + margin
    npts = np.ceil((hi - lo) / spacing_bohr).astype(int) + 1
    return lo, spacing_bohr, npts


def evaluate(mol, dm, lo, delta, npts, chunk=None, mem_mb=400, label=""):
    """Wertet rho(r) und V(r) auf dem Gitter aus.

    V(r) wird ueber ``int1e_grids`` blockweise berechnet - das ist die
    PySCF-Routine, die genau fuer Potentialauswertungen an vielen Punkten
    gedacht ist. Blockweise, weil die Zwischenmatrix sonst der Speicherfresser
    waere: pro Gitterpunkt eine volle Basisfunktionsmatrix.
    """
    from pyscf.dft import numint

    # Blockgroesse aus dem Speicherbedarf ableiten, nicht raten.
    # int1e_grids liefert ein Feld der Form (punkte, nao, nao) - bei 193
    # Basisfunktionen sind das schon 0.28 MB PRO Gitterpunkt. Mit einem festen
    # Block von 20000 Punkten wollte diese Funktion 6 GB anfordern und wurde
    # vom Betriebssystem abgeschossen. Deshalb: Block so waehlen, dass die
    # Zwischenmatrix unter mem_mb bleibt.
    nao = mol.nao
    if chunk is None:
        chunk = max(64, int(mem_mb * 1024**2 / (nao * nao * 8)))

    ax = [lo[i] + np.arange(npts[i]) * delta for i in range(3)]
    # Reihenfolge wie Turbomole: x laeuft am schnellsten, z am langsamsten
    Z, Y, X = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")
    grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)
    total = grid.shape[0]

    rho = np.empty(total)
    esp = np.empty(total)

    charges = mol.atom_charges().astype(float)
    nuclei = mol.atom_coords()                       # Bohr

    if label:
        print(f"    {nao} Basisfunktionen -> Block von {chunk:,} Punkten "
              f"({chunk * nao * nao * 8 / 1024**2:.0f} MB Zwischenmatrix)")

    t0 = time.time()
    for start in range(0, total, chunk):
        pts = grid[start:start + chunk]

        # Elektronendichte
        ao = numint.eval_ao(mol, pts)
        rho[start:start + chunk] = numint.eval_rho(mol, ao, dm)

        # Elektrostatisches Potential: Kerne minus Elektronen
        d = np.linalg.norm(pts[:, None, :] - nuclei[None, :, :], axis=2)
        d[d < 1e-12] = 1e-12
        v_nuc = (charges[None, :] / d).sum(axis=1)

        # int1e_grids liefert <i|1/|r-C||j> mit POSITIVEM Vorzeichen; der
        # elektronische Beitrag zum Potential muss also abgezogen werden.
        # Kontrolle: weit ausserhalb eines neutralen Molekuels muss V gegen 0
        # gehen. Bei HCl in 40 Bohr Abstand: v_nuc = +0.4516, v_ele = +0.4512,
        # Differenz +0.0004 - Summe waere +0.9028 und damit offensichtlich falsch.
        ints = mol.intor("int1e_grids", grids=pts)
        v_ele = np.einsum("pij,ij->p", ints, dm)

        esp[start:start + chunk] = v_nuc - v_ele

        done = min(start + chunk, total)
        if label:
            pct = 100.0 * done / total
            eta = (time.time() - t0) / done * (total - done)
            sys.stdout.write(f"\r    {label}: {pct:5.1f} %  "
                             f"(noch ~{eta / 60:.1f} min)   ")
            sys.stdout.flush()
    if label:
        sys.stdout.write(f"\r    {label}: 100.0 %  "
                         f"({total:,} Punkte in {(time.time()-t0)/60:.1f} min)\n")

    return rho, esp


# ----------------------------------------------------------------------------
# Turbomole-pointval-Format schreiben
# ----------------------------------------------------------------------------

def write_pointval(path, values, lo, delta, npts, quantity, title):
    """Schreibt ein Gitter im Turbomole-``pointval``-Format.

    Bewusst dasselbe Format wie die echten Daten, inklusive der vollen
    Koordinaten pro Zeile und der Reihenfolge x-schnellst. Nur so durchlaeuft
    der Testfall auch xyzToCube.py.
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
# Ein Molekuel komplett
# ----------------------------------------------------------------------------

def make_case(name, smiles, root, spacing, margin, basis, method):
    print(f"\n[{name}]  {smiles}")
    outdir = os.path.join(root, name)

    atoms, molfile = build_geometry(smiles, name, outdir)
    print(f"    Geometrie: {len(atoms)} Atome -> {os.path.basename(molfile)}")

    t0 = time.time()
    mol, dm, energy = run_scf(atoms, basis=basis, method=method)
    print(f"    SCF ({method}/{basis}): E = {energy:.6f} Hartree, "
          f"{mol.nao} Basisfunktionen, {time.time()-t0:.1f} s")

    lo, delta, npts = make_grid(atoms, spacing_bohr=spacing, margin_ang=margin)
    total = int(np.prod(npts))
    print(f"    Gitter: {npts[0]} x {npts[1]} x {npts[2]} = {total:,} Punkte, "
          f"delta = {delta} Bohr")

    rho, esp = evaluate(mol, dm, lo, delta, npts, label="rho und V")

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
        description="Testdatensaetze (td.xyz/tp.xyz) aus SMILES erzeugen.")
    p.add_argument("--smiles", help="SMILES eines einzelnen Molekuels")
    p.add_argument("--name", help="Ordner- und Dateiname dazu")
    p.add_argument("--preset", action="store_true",
                   help="die beiden eingebauten Testfaelle erzeugen")
    p.add_argument("--outdir", default="../sandbox",
                   help="Wurzelverzeichnis fuer die Molekuelordner")
    p.add_argument("--spacing", type=float, default=0.25,
                   help="Gitterabstand in Bohr (Standard 0.25)")
    p.add_argument("--margin", type=float, default=3.5,
                   help="Rand um das Molekuel in Angstrom (Standard 3.5)")
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--method", default="hf", choices=["hf", "b3lyp"])
    args = p.parse_args(argv)

    print("=" * 70)
    print("CreateTpTdFromSmiles.py - Testdaten fuer den ESP-Workflow")
    print("=" * 70)
    print("Hinweis: Kraftfeld-Geometrie, eigenes Rechenniveau.")
    print("Das sind Testfixtures - keine Referenzdaten fuer Vergleichstabellen.")

    cases = list(PRESETS) if args.preset else []
    if args.smiles:
        if not args.name:
            raise SystemExit("--smiles braucht auch --name")
        cases.append((args.name, args.smiles))
    if not cases:
        raise SystemExit("Nichts zu tun: --preset oder --smiles/--name angeben.")

    made = []
    for name, smiles in cases:
        made.append(make_case(name, smiles, args.outdir, args.spacing,
                              args.margin, args.basis, args.method))

    print("\nFertig. Erzeugte Ordner:")
    for m in made:
        print(f"  {m}")
    print("\nWeiter mit:")
    print("  cd ../scripts && python run_all.py --root ../sandbox --stride 1")
    return 0


if __name__ != "CreateTpTdFromSmiles":
    _argv = sys.argv[1:]
    for _i, _a in enumerate(_argv):
        if _a.endswith("CreateTpTdFromSmiles.py"):
            _argv = _argv[_i + 1:]
            break
    main(_argv)
