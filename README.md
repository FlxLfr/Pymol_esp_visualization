# ESP Visualization

**Standard Operating Procedure — How to Visualize Molecular Electrostatic Potential Data**

A reproducible, scriptable workflow that turns Turbomole `pointval` output into
publication-quality images of the molecular electrostatic potential (ESP) mapped
onto an electron-density isosurface.

| | |
|---|---|
| Input | Turbomole `pointval` grids (`td.xyz`, `tp.xyz`) + a structure file |
| Output | Gaussian cube files, a standard set of PNG images, a CSV of surface ESP statistics |
| Software | Python 3 + NumPy + PyMOL (open source), all free |
| Manual steps | none — the whole pipeline is two commands |

<p align="center">
  <img src="reference/4-bromacetophenon/images/4-bromacetophenon_pi.png" width="32%" alt="pi face">
  <img src="reference/4-bromacetophenon/images/4-bromacetophenon_sigma.png" width="32%" alt="sigma hole">
  <img src="reference/4-bromacetophenon/images/4-bromacetophenon_edge.png" width="32%" alt="in-plane profile">
</p>
<p align="center">
  <img src="reference/4-bromacetophenon/images/4-bromacetophenon_colorbar.png" width="42%" alt="colour scale">
</p>

<p align="center"><em>4-Bromoacetophenone. Left: π face, with the carbonyl
oxygen as the deep red region at the top. Centre: view along the C–Br axis — the
blue σ-hole cap inside the red belt. Right: in-plane profile. One molecule
showing both features the workflow has to get right: a strongly negative
functional group and an anisotropic halogen.</em></p>

---

## Contents

1. [Why these images look the way they do](#1-why-these-images-look-the-way-they-do)
2. [Installation](#2-installation)
3. [Input files and formats](#3-input-files-and-formats)
4. [Step 1 — Convert the grids to cube](#4-step-1--convert-the-grids-to-cube)
5. [Step 2 — Look at it interactively](#5-step-2--look-at-it-interactively)
6. [Step 3 — Render the standard image set](#6-step-3--render-the-standard-image-set)
7. [Step 4 — Several molecules at once](#7-step-4--several-molecules-at-once)
8. [Choosing the colour scale](#8-choosing-the-colour-scale)
9. [Parameter study](#9-parameter-study)
10. [Results](#10-results)
11. [Repository layout](#11-repository-layout)
12. [Troubleshooting](#12-troubleshooting)
13. [References](#13-references)

---

## 1. Why these images look the way they do

Two conventions are worth stating up front, because they are what make the
images comparable between molecules and defensible in a report.

**The surface is the ρ = 0.001 a.u. electron-density isosurface**, not a
van-der-Waals surface built from tabulated radii. This is the convention
introduced by Politzer, Murray and co-workers: it encloses roughly 96–97 % of
the electronic charge, it responds to the actual electronic structure of the
particular molecule, and it is defined by a single number rather than by a table
of element radii that different programs disagree about.

**The colour is the ESP on that surface, in atomic units (Hartree/e).** Red is
negative (electron-rich, attracts electrophiles), blue is positive. The scale is
symmetric around zero so that white always means "neutral".

For a bromoarene this pair of conventions makes one specific feature visible:
the **σ-hole**. The potential on the bromine is not isotropic — it is *negative*
in a belt perpendicular to the C–Br bond and *positive* in a cap on the extension
of that bond. That positive cap is what allows a halogen bond, and no single
partial charge on the bromine can represent it. The `_sigma` view exists
specifically to show it.

Measured values for the included example, 4-bromoacetophenone:

| quantity | on | a.u. | kcal/(mol·e) | kJ/(mol·e) |
|---|---|---|---|---|
| V<sub>S,min</sub> | carbonyl O | −0.0653 | −41.0 | −171.5 |
| V<sub>S,max</sub> | a ring H | +0.0476 | +29.9 | +125.0 |
| σ-hole | Br, on the C–Br axis | +0.0230 | +14.5 | +60.5 |
| halogen belt | Br, perpendicular | −0.0167 | −10.5 | −43.9 |

Note that all four are different places on the same surface, and that neither
global extremum sits on the bromine — see
[§6](#which-number-describes-the-σ-hole).

---

## 2. Installation

### 2.1 Prerequisite: a working conda

This workflow needs a **functioning conda installation**. That sounds trivial
and is the single most likely thing to cost you an afternoon, so it is spelled
out here.

If you have no conda yet, install
**[Miniforge](https://conda-forge.org/download/)** — it is free, defaults to the
conda-forge channel that this workflow uses, and is small. During installation:

| Option | Setting |
|---|---|
| Create shortcuts | **yes** — this gives you the "Miniforge Prompt" in the Start menu |
| Add to PATH | **no** — the installer is right, it causes conflicts |
| Register as default Python | **no** — leave any existing Python installation alone |
| Clear package cache | yes |

> **Do not rely on a conda that came bundled with another program.** Several
> chemistry packages (Schrödinger PyMOL among them) ship their own conda inside
> their installation folder. It may be incomplete or broken, it is not intended
> to be used as a general package manager, and — worse — VS Code will happily
> discover it and remember it as "the" conda of the system. You then get
> `conda.exe is not a valid application for this operating system platform` on
> every activation, no matter which interpreter you select. See
> [Troubleshooting](#12-troubleshooting) for how to get out of that.

Verify that conda works *before* going further:

```bash
conda --version
```

On Windows, use the **Miniforge Prompt** from the Start menu. To make conda
available in PowerShell as well, run once:

```bash
conda init powershell
```

then close and reopen every shell.

### 2.2 Create the environment

Keeping it separate from `base` means it can be recreated exactly and nothing
else on the machine is disturbed.

```bash
conda env create -f environment.yml
conda activate esp
```

Or explicitly:

```bash
conda create -n esp -c conda-forge python=3.12 pymol-open-source numpy matplotlib
conda activate esp
```

### 2.3 Verify

```bash
python -c "import pymol, numpy, matplotlib; print('ok')"
python -c "import sys; print(sys.executable)"
```

The second line must point *inside* the `esp` environment. If it points at a
system Python instead, the environment was not activated — see
[Troubleshooting](#12-troubleshooting).

---

## 3. Input files and formats

Per molecule, in one folder:

| File | Content | Units |
|---|---|---|
| `td.xyz` | Turbomole `pointval` **total density** grid | Bohr |
| `tp.xyz` | Turbomole `pointval` **total potential** (ESP) grid | Bohr |
| `*.mol` / `*.sdf` / `*.xyz` | molecular structure | Å (default) |

**`td.xyz` and `tp.xyz` are not structure files** despite the extension. They are
ASCII point clouds — one line per grid point, carrying the full coordinates plus
the value:

```
#grid1  start  -15.000000  delta    0.120000  points    251
#electrostatic potential
# cartesian coordinates x,y,z and f(x,y,z)
      -15.00000000   -15.00000000   -15.00000000   -0.00054019
```

At 251³ points that is about 1.25 GB per file. The same information as a
Gaussian cube is roughly 200 MB, because the cube format stores the grid
implicitly.

### Accepted structure formats

Both `xyzToCube.py` and `render_esp.py` accept the same three formats:

| Format | Notes |
|---|---|
| `.xyz` | `Symbol x y z`. With or without the leading atom-count and comment lines — a bare coordinate list is accepted. Unit set by `--struct-unit` (default Å). |
| `.mol` | MDL molfile, V2000 and V3000. **Coordinates come *before* the element symbol**, the reverse of xyz. Always Å. |
| `.sdf` | SD-file; only the first record (up to `$$$$`) is read. Always Å. |

`--struct-unit` applies to `.xyz` only. Molfile coordinates are in Ångström by
definition, so the option is ignored for them.

**Prefer `.mol`/`.sdf` over `.xyz`**: they carry bond information, so PyMOL draws
proper sticks instead of guessing connectivity from distances.

**Give both scripts the same structure file.** `xyzToCube.py` writes the atom
positions into the cube header, `render_esp.py` hands the structure to PyMOL for
the stick model. If you feed them two different files and those files ever
disagree — a different conformer, a different atom order — the skeleton will sit
offset from its own surface, and **nothing will warn you**. One file for both
removes the failure mode entirely.

---

## 4. Step 1 — Convert the grids to cube

```bash
cd scripts
python xyzToCube.py --struct ../path/to/molecule.mol ../path/to/td.xyz ../path/to/tp.xyz --pymol
```

`--struct` takes `.xyz`, `.mol` or `.sdf` — see
[§4](#accepted-structure-formats).

This writes `td.cube`, `tp.cube` and (with `--pymol`) a ready-to-use `esp.pml`
next to the input.

Useful options. The first group changes the cube files, the second only the
generated scene — `--help` shows them separated for that reason:

| Option changing the **cube files** | Effect |
|---|---|
| `--stride 2` | keep every 2nd grid point in each direction → **8× smaller** files |
| `--struct-unit bohr` | structure file is already in Bohr (ignored for `.mol`/`.sdf`) |
| `--outdir DIR` | write the cubes somewhere else |

| Option changing only **`esp.pml`** | Default | Effect |
|---|---|---|
| `--pymol` | off | write the scene at all — without it the four below do nothing |
| `--esp-range` | `auto` | colour range; `auto` derives it from the ESP on the isosurface, exactly as `render_esp.py` does |
| `--pml-iso` | `0.001` | isovalue **drawn in the scene**. Named apart from `render_esp.py`'s `--iso` on purpose: that one moves the measured numbers, this one only the picture |
| `--transparency` | `0.15` | `0` = opaque surface |
| `--rainbow` | off | rainbow ramp in the scene — see [§8.1](#81-the-rainbow-ramp) |

**On `--stride`.** For bromobenzene, decimating the 251³ grid to 126³ leaves
V<sub>S,min</sub> unchanged and shifts V<sub>S,max</sub> by 0.9 % (+19.58 vs.
+19.80 kcal/(mol·e)). The rendered images are indistinguishable, the files shrink
from 201 MB to 26 MB and PyMOL becomes noticeably faster.

**The σ-hole is more sensitive.** It is evaluated by ray casting with
interpolation rather than read off grid points, which makes it far more robust,
but a coarse grid still smooths the isosurface and biases the value low by a few
percent (see [§6](#which-number-describes-the-σ-hole)). So: `--stride 2` while
you are exploring and for the images, full resolution whenever a σ-hole value
goes into a table.

What the converter takes care of, which is where hand-rolled conversions usually
go wrong:

* **Index order.** Turbomole varies *x* fastest, the cube format varies *z*
  fastest. Without reordering you get a transposed, mirrored molecule.
* **Units.** The grid is in Bohr, the structure file is normally in Å (factor
  1.8897). Get this wrong and the molecule floats outside its own surface.

---

## 5. Step 2 — Look at it interactively

Before rendering anything, check that structure and grids actually line up:

```bash
pymol esp.pml
```

or, inside a running PyMOL:

```
cd /path/to/molecule
@esp.pml
```

The script loads the structure and both cubes, builds the ρ = 0.001 isosurface,
maps the ESP onto it and shows the colour ramp. Rotate it. The molecular skeleton
should sit inside its surface, not next to it.

Handy while exploring:

```
turn x, 90          # tip into the ring plane
set transparency, 0     # opaque, strongest colours
set transparency, 0.15  # skeleton shows through (default)
disable espramp     # hide the colour bar
```

---

## 6. Step 3 — Render the standard image set

```bash
cd scripts
python render_esp.py --density ../path/td.cube --esp ../path/tp.cube \
                     --struct ../path/molecule.mol --prefix molecule
```

With no arguments it looks for `td.cube`, `tp.cube` and a structure file in the
current directory. Output goes to `images/`.

Three views are produced, oriented **from the molecular geometry**, not from
PyMOL's `orient`:

| File | View | Shows |
|---|---|---|
| `*_pi.png` | perpendicular to the molecular plane | π system, ring hydrogens |
| `*_sigma.png` | along the C–halogen axis, from outside | **the σ-hole, head on** |
| `*_edge.png` | in the molecular plane | overall profile |

The plane normal comes from the principal axes of the heavy atoms; the second
axis is the carbon–halogen bond. Every molecule therefore lands in the same
orientation automatically — that is what makes the set comparable, and it is why
no manual rotation is needed or wanted. If the molecule carries more than one
halogen, the axis is the C–X bond of the one with the *strongest* σ-hole (see
[Molecules with more than one halogen](#molecules-with-more-than-one-halogen)).

**Two versions of that axis exist, and the difference matters.** The π and edge
views need an axis that lies exactly in the fitted plane, so the three views form
a clean orthogonal triple; the raw C–X vector is projected into the plane for
them. The σ view must not use that projection — it has to follow the *actual*
bond. On a planar molecule the two coincide (4-bromoacetophenone: 0.001°) and
the distinction is invisible. On triazolam the C–Cl21 bond sticks 42.9° out of
the fitted plane, and the σ view was aiming 42.9° past the σ-hole: the picture
showed the molecule roughly edge-on with the chlorine off to the side. Since the
fix, `molecular_frame()` returns both vectors and the σ view uses the unprojected
one. When the tilt exceeds 15°, the console says so.

**For molecules without a halogen** there is no C–X axis, and the script falls
back to the longest principal axis of the heavy atoms. The three views stay a
proper orthogonal triple, but `*_sigma.png` then looks down the long axis of the
molecule rather than at a σ-hole — for paracetamol, for instance, from the
phenol end through the ring towards the acetamido group. The file name is kept
for consistency across an image set; read it as "axial view" in that case.

Alongside the images:

* `*_colorbar.png` — the colour scale as a separate figure (needs matplotlib)
* `*_settings.txt` — every parameter used, including the surface ESP values.
  Keep this next to the figures; it is the record of how they were made.

### Which number describes the σ-hole

**Not `V_S,max`.** On an aryl halide the global maximum of the surface potential
sits on the *ring hydrogens*, not on the halogen. Compare bromobenzene and
iodobenzene:

All values below from the full 251³ grid:

| | chlorobenzene | bromobenzene | iodobenzene |
|---|---|---|---|
| **σ-hole** (cap on the C–X axis) | **+0.0078** a.u. (+4.9 kcal/mol) | **+0.0162** a.u. (+10.2) | **+0.0255** a.u. (+16.0) |
| ring hydrogens = global V<sub>S,max</sub> | +0.0313 a.u. (+19.6) | +0.0315 a.u. (+19.8) | +0.0317 a.u. (+19.9) |
| halogen belt = global V<sub>S,min</sub> | −0.0190 a.u. (−11.9) | −0.0188 a.u. (−11.8) | −0.0169 a.u. (−10.6) |

The three global maxima agree to within 0.4 % — because all three are the same
aromatic C–H. Use them to compare halogen-bond donors and you conclude that
chlorine, bromine and iodine have equally strong σ-holes, which is wrong: the
σ-holes span a factor of 3.3, in the expected order Cl < Br < I.

`render_esp.py` therefore reports both. It names the atom the global extremum
belongs to, and gives the local values in the halogen regions:

```
    V_S,max = +0.0312 a.u.  =  +19.6 kcal/(mol*e)   auf H5
  Lokal am Halogen (Br):
    sigma-Loch  = +0.0126 a.u.  =   +7.9 kcal/(mol*e)   [144 Punkte]
    Guertel     = -0.0188 a.u.  =  -11.8 kcal/(mol*e)   [836 Punkte]
    ! V_S,max liegt auf H5, nicht auf dem Halogen
```

Both values also go into `*_settings.txt` and into the `sigma_hole_au` /
`belt_min_au` columns of `summary.csv`.

In the console the molecule headers are green and halogen symbols are cyan, so
the relevant lines stand out in a long batch run. Colours switch off
automatically when the output is redirected to a file or piped, so log files stay
clean.

To turn them off explicitly, use the flag — it works on both scripts:

```bash
python run_all.py --root ../sandbox --no-color
```

The `NO_COLOR` / `FORCE_COLOR` environment variables are honoured as well
(see [no-color.org](https://no-color.org)). Note that these are *environment
variables*, not arguments, so they are set before the command rather than
appended to it:

```powershell
# PowerShell
$env:NO_COLOR = 1
python run_all.py --root ../sandbox
Remove-Item Env:NO_COLOR
```

```bash
# bash
NO_COLOR=1 python run_all.py --root ../sandbox
```

**The σ-hole is not read off grid points.** It is a peak *on* the C–X axis,
and whether a grid point happens to sit both on that axis and inside the thin
ρ = 0.001 shell is luck. On a decimated bromobenzene grid the best cap point
was 1.14 Bohr off the axis, which measured the flank instead of the summit.

The script therefore casts rays from the halogen into a cone around the axis,
locates where ρ crosses the isovalue along each ray by interpolation, and
evaluates V there — both trilinearly interpolated. That removes the dependence
on grid alignment:

| grid | spacing | point-based | **ray-based (used)** |
|---|---|---|---|
| 42³ | 0.72 Bohr | +2.1 kcal/mol | **+9.4** |
| 126³ (`--stride 2`) | 0.24 Bohr | +7.9 | **+10.1** |
| 251³ (full) | 0.12 Bohr | +11.0 | **+10.2** |

Measured on bromobenzene. The point-based value spans a factor of five across
these grids; the ray-based one varies by 8 %. The two also differ on the *same*
grid, and not by accident: the point-based method takes the best grid point
within a shell of |ρ − 0.001| < 0.00012, and the inner edge of that shell sits
closer to the nucleus where V is higher. For bromobenzene, ρ = 0.00112 gives
+10.9 kcal/mol and ρ = 0.00088 gives +9.3 — the ray method evaluates at exactly
0.00100 and returns +10.2. A residual underestimate remains on very coarse grids because
the interpolated density itself smooths the isosurface, so the script warns
whenever the spacing exceeds 0.30 Bohr. `summary.csv` records which method
produced each value in the `sigma_method` column.

### Molecules with more than one halogen

Two halogens in the same molecule are usually not equivalent. Triazolam is the
obvious example: one chlorine sits on the fused benzodiazepine ring, the other
on the pendant phenyl, and they see completely different electronic
environments. An earlier version of the script simply took the first halogen in
the atom list — which halogen that was depended on nothing but the ordering in
the structure file, and the second one was silently dropped.

Every halogen with a bonded carbon is now evaluated separately: its own C–X
axis, its own ray cone, its own belt. The console prints one block per halogen,
sorted by σ-hole strength, and marks the one the σ view looks along:

```
  Lokal an Br1:   <- Achse der sigma-Ansicht
    sigma-Loch  = +0.0086 a.u.  =   +5.4 kcal/(mol*e)   [interpoliert, 2.9 Grad zur Achse]
    Guertel     = -0.0362 a.u.  =  -22.7 kcal/(mol*e)   [224 Punkte]
  Lokal an Cl4:
    sigma-Loch  = -0.0040 a.u.  =   -2.5 kcal/(mol*e)   [interpoliert, 2.2 Grad zur Achse]
    Guertel     = -0.0362 a.u.  =  -22.7 kcal/(mol*e)   [157 Punkte]
```

The labels are `symbol` + 1-based atom number, so they can be matched against
the structure file directly. `*_settings.txt` gets one line per halogen and
records which one defined the view. In `summary.csv`, `sigma_hole_au` remains
the *strongest* σ-hole so the column stays sortable across a batch;
`sigma_hole_on` names the atom it belongs to, `n_halogens` how many there are,
and `sigma_holes_all` lists every one of them as `label:value` pairs, e.g.
`Br1:0.00862;Cl4:-0.00398`. Nothing is lost. The summary table at the end of a
batch run appends `Br1 of 2` to make clear the printed value is one of several.

Sanity checks: on *p*-dichlorobenzene, where the two chlorines are related by
symmetry, both come out at +0.0074 a.u. On 1-bromo-2-chloroethane they differ as
chemistry demands — Br +5.4, Cl −2.5 kcal/mol, the alkyl chlorine having no
σ-hole worth the name. Molecules with a single halogen are unaffected: the
values and the console output are identical to before.

One caveat for small molecules: the belt region is defined relative to its own
halogen, and if two halogens sit within a few bonds of each other the belts
overlap, so both report the same minimum — visible above, where both belt values
are −0.0362. On molecules where the halogens are far apart, as in triazolam, this
does not arise.

#### The cone must not reach the rest of the molecule

Triazolam exposed a second problem, and this one was not specific to having two
halogens — it had simply never been triggered by the flat test molecules. The
ray method walked outward from the halogen and took the **outermost** crossing
of ρ = 0.001. Triazolam is folded: the cone around the C–Cl21 axis points at the
methyl group on the triazole ring, roughly 4 Å away. The rays left the chlorine
surface, crossed vacuum, and struck the methyl — and its surface was reported as
the σ-hole:

| | before | after |
|---|---|---|
| Cl21 σ-hole | +0.0299 a.u. (**+18.8** kcal/mol), 29.0° off axis | +0.0171 a.u. (**+10.7**), 3.8° |
| Cl21 belt | −0.0670 a.u. (−42.0) — actually the triazole N3/N4 lone pairs | −0.0184 (−11.5) |

+18.8 kcal/mol for an aryl chloride was the tell: chlorobenzene gives +4.9, and
no substituent pattern triples that. The two symptoms in the printed output are
the **off-axis angle** (29° is nearly the cone edge — a real σ-hole peaks within
a few degrees of the axis) and the cap radius, which came out at 3.94 Å instead
of the ~2.0 Å a chlorine surface sits at.

Two changes fix it:

* the ray takes the **first** downward crossing of ρ = 0.001, not the last. The
  ray starts deep inside the halogen's own density, so the first place it drops
  through the isovalue is that halogen's own surface, whatever lies beyond.
* both the ray search and the grid-point cap are cut off at 1.6 × the Bondi van
  der Waals radius of the halogen. The ρ = 0.001 surface sits at roughly
  1.1–1.2 vdW radii, so this leaves plenty of room and excludes everything else.

Single-halogen molecules are unaffected — for 4-bromoacetophenone the σ-hole
moves from +0.02304 to +0.02303 a.u. and the belt is unchanged, because on a
convex molecule the first and the last crossing are the same point. The rendered
images are bit-identical; only the numbers changed.

The lesson is worth keeping: **three flat halobenzenes cannot test a method that
assumes the space beyond the halogen is empty.** A folded drug molecule can. The
same holds for the σ view: the projection of the C–X axis into the fitted plane
was harmless on every planar test molecule and wrong on the first non-planar one.

#### The colour scale can hide a correct σ-hole

Worth knowing before concluding an image is broken. With `--esp-range auto` the
scale is set by the *global* extremes — for triazolam the triazole nitrogens at
−0.084 a.u. On that scale a σ-hole of +0.017 is 20 % of full blue, i.e. almost
white, even though the view is aimed correctly. Rendering the same molecule with
`--esp-range 0.035` makes the cap and its surrounding belt plainly visible. The
number in the console is the evidence for a σ-hole; the picture at the automatic
scale is not necessarily.

`--stride 2` is fine for the images and for V<sub>S,min</sub> /
V<sub>S,max</sub> — those change by about 1 %. For σ-hole values that go into a
table, use the full grid.

Options worth knowing:

| Option | Default | Effect |
|---|---|---|
| `--esp-range` | `auto` | `auto` or a fixed value in a.u. — see [§8](#8-choosing-the-colour-scale) |
| `--transparency` | `0.15` | `0` = opaque, strongest colours; above ~0.3 the profile views become unreadable |
| `--backgrounds white black` | `white` | render each view on both backgrounds |
| `--views pi sigma` | all three | subset of views |
| `--width / --height / --dpi` | 2000 / 1600 / 300 | image size |
| `--iso` | `0.001` | density isovalue |
| `--buffer` | `2.4` | margin around the molecule, Å |
| `--rainbow` | off | rainbow ramp instead of red–white–blue — see [§8.1](#81-the-rainbow-ramp) |

**Do not use the PyMOL launcher unless you have to.** `python render_esp.py …`
loads no `pymolrc`, so nobody's personal start-up file can silently change a
setting. If you do use the launcher, both the `--` separator and `-k` are
required:

```bash
pymol -ckq render_esp.py -- --prefix molecule
```

---

## 7. Step 4 — Several molecules at once

Put each molecule in its own folder under a common root — for your own data that
is `sandbox/`, which git ignores:

```
sandbox/
├── bromobenzene/   td.xyz  tp.xyz  bromobenzene.mol
├── iodobenzene/    td.xyz  tp.xyz  iodobenzene.mol
└── chlorobenzene/  td.xyz  tp.xyz  chlorobenzene.mol
```

Then:

```bash
cd scripts
python run_all.py --root ../sandbox --two-pass
```

Conversion runs at **full grid resolution** by default, matching
`xyzToCube.py`. Add `--stride 2` for a faster first pass — eight times smaller
cubes, visually identical images, still ~140 grid points in the σ-hole cap:

```bash
python run_all.py --root ../sandbox --two-pass --stride 2
```

`--stride` only matters while cube files are being *created*. If `td.cube` and
`tp.cube` already exist they are reused unchanged, and the flag does nothing —
use `--force-convert` to rebuild them.

To pick individual molecules out of a larger root, use `--only`; simple
wildcards work:

```bash
python run_all.py --root ../sandbox --only paracetamol chlormethan --two-pass
python run_all.py --root ../sandbox --only "*benzol"
```

**Do not put a common colour scale across data of different provenance.**
`--two-pass` gives every molecule in the run the same scale, which is only
meaningful if they were computed the same way — same geometry optimisation,
same method, same basis set. Mixing the provided Turbomole data with grids
generated by `tools/CreateTpTdFromSmiles.py` in one `--two-pass` run would
imply a comparability that does not exist. Use `--only` (or separate root
folders) to keep the groups apart, and run `--two-pass` within each group.

There is a practical side to this as well: the halobenzenes need ±0.035 a.u.,
paracetamol ±0.090. On a shared scale the halobenzenes would come out almost
colourless — a σ-hole of 0.0070 is 8 % of a ±0.090 range.

Called **without arguments**, `run_all.py` runs on `reference/` instead. That is
the smoke test: it exercises the whole pipeline on data that is known to work, so
you can tell an installation problem from a data problem before touching your own
files.

```bash
python run_all.py
```

It converts `td.xyz`/`tp.xyz` to cube files, renders, and writes to
`reference/*/images_check/` and `reference/summary_check.csv` — never to the
committed `images/`, so you can compare side by side. Everything it produces is
git-ignored. Your run should reproduce, for 4-bromoacetophenone:

| | expected |
|---|---|
| V<sub>S,min</sub> | −0.0638 a.u. on O3 |
| V<sub>S,max</sub> | +0.0469 a.u. on H14 |
| σ-hole | +0.0221 a.u. |
| colour range | ±0.065 a.u. |

The images will look coarser than the committed ones, and the numbers differ
slightly from those in [§1](#1-why-these-images-look-the-way-they-do): the smoke
test runs on the decimated 0.75 Bohr grid, the reference images were rendered at
0.25 Bohr. The script says so itself, warning that 0.75 Bohr is too coarse for a
σ-hole value.

This converts what needs converting, renders every molecule, writes an
`esp.pml` next to each molecule's cube files so you can open the scene
interactively, and collects `summary.csv` with V<sub>S,min</sub>,
V<sub>S,max</sub> and the σ-hole for each — in a.u., kcal/(mol·e) and
kJ/(mol·e).

The generated `esp.pml` always carries the colour scale that was actually used
for that molecule's images. After a `--two-pass` run it holds the common scale,
so what you see interactively matches the figure set.

`--two-pass` first renders every molecule on its own automatic scale, then
re-renders all of them using the largest range it saw, so the final set is
directly comparable. That is the recommended mode for a figure set.

---

## 8. Choosing the colour scale

This is the one decision that cannot be automated away, because it depends on
what the figure is for.

**One molecule, maximum detail** → `--esp-range auto`. The range is derived from
the ESP actually present on the ρ = 0.001 shell and rounded up to a clean value.

**Several molecules to be compared** → one fixed range for all of them.
A red patch in figure A must mean the same potential as a red patch in figure B,
otherwise the comparison is misleading. `--two-pass` does this for you; manually,
run once with `auto`, take the largest reported value and re-run everything with
`--esp-range <value>`.

Whatever you choose, **state the range in the figure caption** and ship
`*_colorbar.png` with the figures. An ESP figure without its scale is
uninterpretable.

For the halobenzenes the automatic range comes out at **±0.035 a.u.**
(±92 kJ/(mol·e)), stable across grid resolutions; 4-bromoacetophenone needs
±0.070 because of the carbonyl oxygen.

### 8.1 The rainbow ramp

`--rainbow` switches the colour ramp from red–white–blue to a rainbow. It works
on all three scripts:

```bash
python run_all.py --root ../sandbox --rainbow
python render_esp.py --density td.cube --esp tp.cube --struct x.mol --rainbow
python xyzToCube.py --struct x.mol td.xyz tp.xyz --pymol --rainbow
```

| ramp | levels | colours |
|---|---|---|
| default | −rng, 0, +rng | red, white, blue |
| `--rainbow` | −rng, −rng/2, 0, +rng/2, +rng | red, yellow, green, cyan, blue |

**The convention is deliberately not inverted.** Red stays negative and blue
stays positive in both ramps; the rainbow only inserts yellow, green and cyan in
between. Some programs run the rainbow the other way round (blue = negative),
which makes the two image sets impossible to place side by side — the point of
having the option at all.

**The two sets do not overwrite each other.** A rainbow run writes its own
files, so the standard set survives:

```
images/<mol>_pi.png          images/<mol>_rainbow_pi.png
images/<mol>_colorbar.png    images/<mol>_rainbow_colorbar.png
images/<mol>_settings.txt    images/<mol>_rainbow_settings.txt
esp.pml                      esp_rainbow.pml
```

`*_settings.txt` records the ramp in a `Farbrampe` line, and `summary.csv` has a
`colormap` column with `redblue` or `rainbow`, so a figure can always be traced
back to the ramp that produced it.

**When to use which.** Red–white–blue spends almost no colour resolution near
zero: everything weakly polar comes out white. The rainbow resolves exactly that
region, which is useful for looking at the π face of an aromatic ring or at a
saturated backbone.

The price is that a rainbow is not perceptually uniform. The eye reads the
yellow–green and green–cyan transitions as edges even where the potential
changes smoothly, and green in the middle looks like a state of its own rather
than "neutral". For quantitative statements — σ-hole, belt, V<sub>S,min</sub>,
V<sub>S,max</sub> — red–white–blue is the more honest picture, and it is also
what the ESP literature uses. Recommended practice: red–white–blue as the main
figure, rainbow alongside where fine structure is the point.

---

## 9. Parameter study

Every default in this pipeline was chosen against a measurement, not by taste.
This section collects those measurements so the choices can be checked, and so
it is clear which parameters are cosmetic and which change the result.

The short version: **the isovalue is the only parameter that changes the
physics.** Grid resolution costs accuracy slowly, the colour range and
transparency change nothing but the picture, and the σ-hole search parameters
matter only through the two failure modes described in
[§6](#6-step-3--render-the-standard-image-set).

### 9.1 Isovalue ρ — the one that decides the answer

4-bromoacetophenone, full 114×86×80 grid, everything else at defaults:

| ρ / a.u. | shell points | V<sub>S,min</sub> | V<sub>S,max</sub> | σ-hole | σ-hole kcal/mol | belt |
|---|---|---|---|---|---|---|
| 0.0005 | 4901 | −0.0596 | +0.0408 | +0.0159 | 10.0 | −0.0151 |
| 0.0008 | 4925 | −0.0640 | +0.0463 | +0.0203 | 12.7 | −0.0162 |
| **0.0010** | **4707** | **−0.0653** | **+0.0476** | **+0.0230** | **14.5** | **−0.0167** |
| 0.0015 | 4743 | −0.0703 | +0.0539 | +0.0293 | 18.4 | −0.0176 |
| 0.0020 | 4666 | −0.0730 | +0.0584 | +0.0354 | 22.2 | −0.0179 |
| 0.0040 | 4490 | −0.0816 | +0.0748 | +0.0568 | 35.7 | −0.0171 |

From 0.0005 to 0.004 the σ-hole grows by a factor of 3.6. That is not noise —
a larger isovalue means a surface closer to the nuclei, where the positive
nuclear contribution has been screened less. Every value on such a surface is
larger in magnitude, on both signs.

The consequence is that **a σ-hole value without its isovalue is meaningless**,
and two values computed at different isovalues cannot be compared at all. This
is why 0.001 a.u. is not a tunable here: it is the Bader/Politzer convention
(see [§1](#1-why-these-images-look-the-way-they-do)) and the only value the rest
of the literature can be read against. `--iso` exists to reproduce someone
else's choice, not to improve on this one.

Note that the shell point count barely moves across the whole range. The count
is therefore no indication that anything changed — a convergence check that
looks at "enough points" would have passed at every one of these settings.

### 9.2 Grid resolution — `--stride`

Same molecule, isovalue fixed at 0.001, cubes rebuilt from the same pointval
files at increasing decimation:

| `--stride` | grid | Δ / Bohr | cubes | V<sub>S,min</sub> | V<sub>S,max</sub> | σ-hole | kcal/mol | belt |
|---|---|---|---|---|---|---|---|---|
| **1** | 114×86×80 | **0.25** | 20.7 MB | −0.0653 | +0.0476 | **+0.0230** | 14.5 | −0.0167 |
| 2 | 57×43×40 | 0.50 | 2.6 MB | −0.0636 | +0.0462 | +0.0224 | 14.1 | −0.0166 |
| 3 | 38×29×27 | 0.75 | 0.8 MB | −0.0638 | +0.0469 | +0.0221 | 13.9 | −0.0159 |
| 4 | 29×22×20 | 1.00 | 0.3 MB | −0.0634 | +0.0428 | +0.0214 | 13.4 | — |

Cost falls by a factor of 70 across this table; the σ-hole falls by 7 %. The
degradation is smooth and one-sided — always low, never high — because the
interpolated density smooths the isosurface. At stride 4 the belt can no longer
be measured at all: fewer than five surface points survive in the belt region,
and the script prints a dash rather than a number.

Two things make this table look better than it should. First, these are
**ray-based** values; the point-based method on the same grids spans a factor of
five (see the table in [§6](#6-step-3--render-the-standard-image-set)). Second,
this molecule is small — for the 251³ bromobenzene grids the full-resolution
cubes are 205 MB and stride 2 is what makes the workflow usable at all.

Practical rule: **stride 2 for looking, stride 1 for numbers that go in a
table.** The script warns above 0.30 Bohr spacing for exactly this reason.

### 9.3 Colour range — visual only

`--esp-range` never touches a computed value; it only maps numbers to colours.
It is in this section because getting it wrong makes a correct calculation look
wrong. Ranges the automatic mode selects:

| molecule | range | set by |
|---|---|---|
| halobenzenes | ±0.035 | the halogen belt |
| 4-bromoacetophenone | ±0.070 | the carbonyl oxygen |
| triazolam | ±0.085 | the triazole nitrogens |
| paracetamol | ±0.090 | the phenol oxygen |

A σ-hole of +0.017 is 49 % of full blue at ±0.035 and 20 % at ±0.085 — visible
in the first case, nearly white in the second, from identical data. Whenever the
molecule carries a group far more polar than the halogen, the automatic range is
set by that group and the halogen region is washed out. See
[§8](#8-choosing-the-colour-scale) for when a common scale is legitimate.

### 9.4 σ-hole search parameters

| parameter | value | why |
|---|---|---|
| cone half-angle | 36.9° (`cone_cos = 0.80`) | wide enough that the maximum is found interior to the cone, narrow enough to exclude the belt. A maximum near the edge is a warning sign, not a result — the console prints the off-axis angle for this reason |
| rays per halogen | 400, Fibonacci spiral | even coverage of the cap without the pole clustering of spherical coordinates |
| step along a ray | 0.02 Bohr | far below the grid spacing; the crossing radius is then refined by linear interpolation anyway |
| radius cut-off | 1.6 × Bondi vdW | the ρ = 0.001 surface sits at 1.1–1.2 vdW radii; beyond that the ray is looking at another part of the molecule |
| belt half-angle | ±69.5° (`belt_cos = 0.35`) | the belt is broad; a narrow band would sample only its rim |
| belt radius | 1.5 × mean cap radius | keeps the belt on its own halogen — with the caveat for adjacent halogens noted in [§6](#6-step-3--render-the-standard-image-set) |

The first two of these were tuned; the rest follow from the geometry. What
actually cost the most work was not any of these numbers but the two structural
mistakes — taking the outermost isosurface crossing, and projecting the C–X axis
into the fitted plane — both of which are documented where they arose.

### 9.5 Rendering

| parameter | value | note |
|---|---|---|
| `transparency` | 0.15 | enough to see the stick model through the surface, little enough that the colour still reads. Purely a visual judgement; no number depends on it |
| `transparency_mode` | 2 | PyMOL's back-face-aware mode; without it the far side of the surface bleeds through and the colours mix |
| `surface_quality` | 1 | the isosurface triangulation, not the data |
| `orthoscopic` | on | no perspective, so two molecules rendered at the same zoom are directly comparable |
| image size | 2000 × 1600, 300 dpi | large enough for a full-width figure in a report at ~17 cm |

### 9.6 Test-data generator

`tools/CreateTpTdFromSmiles.py` has its own parameters, documented with their
reasoning in `tools/README.txt`. The one measurement worth repeating here is the
cost scaling, because it decides what is practical: triazolam (35 atoms, 390
basis functions) at the default 0.25 Bohr and 3.5 Å margin needs 1.41 million
grid points and about 26 minutes for the grid evaluation, against roughly 4
minutes for 4-bromoacetophenone. The two factors multiply — 1.8× more points and
4× more work per point. `--spacing 0.30 --margin 2.5` brings that back to about
10 minutes, at the price of landing at the grid-spacing warning threshold.

---

## 10. Results

`results/` holds the image sets that this project actually delivers — seven
molecules, each with the three standard views, its colour bar and its
`*_settings.txt`. Unlike `sandbox/`, this folder is committed: the images are the
deliverable, and `*_settings.txt` next to each one is the record of how it was
made.

```
results/
├── chlorbenzol/          chlorobenzene       ┐
├── brombenzol/           bromobenzene        │ provided Turbomole data
├── iodbenzol/            iodobenzene         ┘
├── chlormethan/          chloromethane       ┐
├── 4-bromacetophenon/    4-bromoacetophenone │ generated with
├── paracetamol/          paracetamol         │ tools/CreateTpTdFromSmiles.py
└── halcion/              triazolam           ┘
```

### 10.1 Surface ESP values

All at ρ = 0.001 a.u.; σ-hole and belt from the ray method. Atom labels are
1-based indices into the respective structure file.

| molecule | grid | Δ/Bohr | range | V<sub>S,min</sub> | V<sub>S,max</sub> | σ-hole | kcal/mol | belt |
|---|---|---|---|---|---|---|---|---|
| chlorobenzene | 251³ | 0.12 | ±0.035 | −0.0190 (Cl12) | +0.0313 (H9) | +0.0078 | +4.9 | −0.0190 |
| bromobenzene | 251³ | 0.12 | ±0.035 | −0.0188 (Br12) | +0.0315 (H7) | +0.0162 | +10.2 | −0.0188 |
| iodobenzene | 251³ | 0.12 | ±0.035 | −0.0169 (I12) | +0.0317 (H5) | +0.0255 | +16.0 | −0.0169 |
| chloromethane | 40×38×38 | 0.40 | ±0.035 | −0.0283 (Cl1) | +0.0323 (H3) | **−0.0083** | **−5.2** | −0.0283 |
| 4-bromoacetophenone | 114×86×80 | 0.25 | ±0.070 | −0.0653 (O3) | +0.0476 (H14) | +0.0230 | +14.5 | −0.0167 |
| triazolam | 132×101×106 | 0.25 | ±0.085 | −0.0843 (N4) | +0.0543 (H27) | +0.0171 (Cl21) | +10.7 | −0.0184 |
| " | | | | | | +0.0144 (Cl11) | +9.0 | −0.0110 |
| paracetamol | 122×67×86 | 0.25 | ±0.090 | −0.0737 (O3) | +0.0898 (H20) | — | — | — |

What the set is meant to show:

* **The halobenzene series** is the trend Cl < Br < I, a factor of 3.3 across
  the three — while their global V<sub>S,max</sub> agree to within 0.4 %,
  because that maximum is the same aromatic C–H in all three. This is the
  central argument of [§6](#which-number-describes-the-σ-hole).
* **Chloromethane** is the negative control: an *alkyl* chloride has no σ-hole
  at all. The value on the C–Cl axis is −5.2 kcal/mol — still negative. The
  σ-hole is a property of the C–X bond's electronic environment, not of the
  halogen alone.
* **4-Bromoacetophenone** separates the belt from the global minimum for the
  first time: V<sub>S,min</sub> sits on the carbonyl oxygen at −41.0, the belt
  on the bromine at −10.5. In the halobenzenes those two lines were always the
  same number.
* **Paracetamol** is the halogen-free case — no σ-hole block, orientation from
  the principal axes, dash in the table.
* **Triazolam** is the two-halogen case, and the two chlorines differ by 19 %
  despite being the same element: +10.7 on the pendant 2′-phenyl chlorine
  against +9.0 on the fused-ring one.

### 10.2 What these numbers may and may not be compared with

**The first three rows and the last four are not comparable with each other.**
Chlorobenzene, bromobenzene and iodobenzene come from the provided Turbomole
calculation. The other four were generated with `tools/CreateTpTdFromSmiles.py`
at HF/def2-SVP on an MMFF94 geometry — a different method, a different basis
set, and a force-field geometry rather than an optimised one. Any of those three
differences moves V<sub>S,min</sub> and V<sub>S,max</sub> by more than the
effects being discussed. See `tools/README.txt`, section "LIMITATIONS".

Within each group the comparison is sound: the three halobenzenes were computed
identically, and the four generated sets share method, basis and grid spacing.

The colour scales tell the same story. All four molecules at ±0.035 can be
compared *by eye*; 4-bromoacetophenone, triazolam and paracetamol each carry
their own range because a shared one would render the halogen regions
colourless. Always read the colour bar shipped with the image, never the colours
alone.

Two further caveats:

* **Chloromethane's grid is 0.40 Bohr**, above the 0.30 Bohr threshold at which
  the script warns. Per §9.2 that biases the value low by a few percent. It does not affect the conclusion — the sign is
  what matters here, and −5.2 is nowhere near zero.
* The maximum density value in the provided data suggests **bromine was treated
  all-electron while iodine used an effective core potential**. That is a
  question for the supervisor; if so, the Br/I step in the table is not a pure
  basis-set-consistent comparison.

### 10.3 Provenance of the images

> **Note.** `chlorbenzol/`, `chlormethan/` and `iodbenzol/` were rendered with an
> earlier version of `render_esp.py` — their `*_settings.txt` still uses the
> single-halogen layout (`sigma-Loch (Cl)` and a separate `Halogenguertel` line)
> rather than the per-halogen one. The **numbers are unaffected**: all three
> molecules carry one halogen and are planar or axially symmetric, so neither the
> outermost-crossing fix nor the true-axis fix changes anything for them, and the
> values above match a current run. Still, a deliverable set should come from one
> script version. Re-rendering them is one command and the cube files are
> unchanged:
>
> ```bash
> python run_all.py --root ../sandbox --only chlorbenzol chlormethan iodbenzol
> ```

---

## 11. Repository layout

```
esp_visualization/
├── README.md                     this document (the SOP)
├── environment.yml               conda environment
├── .gitignore
├── scripts/
│   ├── xyzToCube.py              Turbomole pointval -> Gaussian cube
│   ├── render_esp.py             standard image set from cube files
│   ├── run_all.py                batch driver + summary.csv
│   ├── ansi.py                   console colours (no dependencies)
│   └── esp.pml                   interactive PyMOL scene
├── reference/                    known-good example — output, not input
│   ├── summary.csv
│   └── 4-bromacetophenon/
│       ├── 4-bromacetophenon.mol
│       ├── td.xyz                raw pointval grids, decimated to 0.75 Bohr
│       ├── tp.xyz
│       └── images/               reference images (rendered at 114×86×80)
├── results/                      the delivered image sets — see §11
│   ├── chlorbenzol/  brombenzol/  iodbenzol/       provided Turbomole data
│   ├── chlormethan/  4-bromacetophenon/
│   ├── paracetamol/  halcion/                      generated test data
│   └── <molecule>/               *_pi.png  *_edge.png  *_sigma.png
│                                 *_colorbar.png  *_settings.txt
├── tools/                        test-data generator (own environment)
├── docs/                         exported PDF of this SOP
└── sandbox/                      your own data and experiments, not tracked
```

**`results/` is committed, `sandbox/` is not.** The images in `results/` are the
deliverable and are small enough to track (PNG + text, a few hundred KB each);
the cube files and pointval grids they were made from stay in `sandbox/` and are
ignored. Each result folder carries its own `*_settings.txt`, so an image never
travels without the parameters that produced it.

**`reference/` and `sandbox/` do different jobs.** `reference/` holds a
known-good example: the images this workflow is supposed to produce, the
parameters that produced them, and a small decimated dataset to reproduce them
with. It is committed, and you do not edit it. `sandbox/` is where your own
molecules and the large raw data live; git ignores it entirely. If a run goes
wrong, `reference/` tells you whether the problem is your installation or your
data.

**Large files are deliberately not tracked.** `.gitignore` excludes `*.cube`
and the Turbomole `td.xyz`/`tp.xyz` grids — a full-resolution cube is 201 MB and
GitHub rejects anything above 100 MB. Regenerate them from the raw data with
`xyzToCube.py`.

**The reference dataset is an exception, and it ships as raw `pointval` files,
not as cubes.** That is deliberate: a smoke test that starts from ready-made
cube files would skip `xyzToCube.py` — the unit conversion and the index
reordering — which is exactly the step most likely to break. Starting from
`td.xyz`/`tp.xyz` exercises the entire chain.

They are decimated to 0.75 Bohr (~2.3 MB each instead of 61 MB, and the original
Turbomole grids are 1.25 GB). The code path does not care about the grid size,
so a coarse grid tests it just as well.

```bash
cd scripts
python run_all.py            # converts, renders, writes summary_check.csv
```

The images in `reference/4-bromacetophenon/images/` were rendered from the full
114×86×80 grid, so they are smoother than what the decimated grid produces, and
the values differ slightly: the reference grid gives V<sub>S,min</sub> = −0.0638
and a σ-hole of +0.0221 against −0.0653 and +0.0230 at full resolution. Use it to
confirm the pipeline runs, not to read numbers off — the script says as much,
warning that 0.75 Bohr is too coarse for a σ-hole value.

Use `sandbox/` for experiments and large data; it is ignored by git, so the
scripts in `scripts/` stay the single source of truth. Do not keep a second copy
of the scripts elsewhere — that is exactly how two versions drift apart.

---

## 12. Troubleshooting

**`ContourSurfVolume: VTKm not available, falling back to internal implementation`**
Harmless, and it appears on every run. VTK-m is an optional parallel contouring
backend that the conda-forge PyMOL build is not compiled with. PyMOL uses its
built-in marching-tetrahedra routine instead; same isosurface, slightly slower.
Nothing to fix.

**The script runs under `pymol -cq` but nothing happens, no error.**
The `--` separator is missing, so PyMOL swallowed the arguments. Use
`pymol -ckq render_esp.py -- --prefix molecule`, or simply
`python render_esp.py --prefix molecule`.

**The molecule floats next to its surface instead of inside it.**
Unit mismatch. The grid is in Bohr, the structure file is probably in Å. Check
the `--struct-unit` setting.

**The molecule looks mirrored or transposed.**
Grid index order. Turbomole varies *x* fastest, cube varies *z* fastest.
`xyzToCube.py` handles this; a hand-written converter usually does not.

**The molecule is tiny in the middle of a large empty image.**
Something zoomed on the isosurface object. A PyMOL isosurface carries the extent
of the *entire grid box*, not of the visible surface. Zoom on the molecule
object instead.

**Colours look washed out, or the profile view is unreadable.**
Transparency too high. Above ~0.3 you are looking through the whole molecule and
seeing the far side's colours mixed in. Use `0.15`, or `0` for maximum contrast.

**No `*_colorbar.png` was produced.**
matplotlib is missing in the active environment: `conda install -c conda-forge
matplotlib`. The molecule images are unaffected — matplotlib only draws the
scale bar.

**`conda: command not found` in `cmd.exe` although it works in VS Code.**
Conda is installed but not on the system PATH; VS Code activates it explicitly.
Run `conda init powershell` (or `conda init cmd.exe`) once from a shell where
conda *does* work — on Windows, the Miniforge Prompt.

**`conda.exe is not a valid application for this operating system platform`.**
The conda being invoked is a broken one bundled with another program. Find out
which one is configured:

```powershell
$env:CONDA_EXE
```

If that prints a path inside some application's folder rather than your
Miniforge installation, that is the culprit. It can come from three places —
check them in this order:

1. **A stale PowerShell profile.** `Test-Path $PROFILE`; if `True`, open it with
   `notepad $PROFILE` and delete the block referring to the foreign conda.
2. **A permanent environment variable.**
   `[Environment]::GetEnvironmentVariable("CONDA_EXE","User")` — clear it with
   `[Environment]::SetEnvironmentVariable("CONDA_EXE",$null,"User")`.
3. **VS Code.** If the variable is set *only* inside the VS Code terminal and
   nowhere else, the Python extension is responsible. Open user settings
   (`Ctrl+Shift+P` → "Preferences: Open User Settings (JSON)") and point it at
   the right conda:

   ```json
   "python.condaPath": "C:\\Users\\<you>\\miniforge3\\Scripts\\conda.exe"
   ```

   Then check `.vscode/settings.json` in the repository for a stale
   `condaPath` or `defaultInterpreterPath`, and reload the window.

Note that selecting the right interpreter is **not** sufficient: the interpreter
decides *which* environment is activated, `CONDA_EXE` decides *what does the
activating*.

**The VS Code terminal uses the wrong Python.**
`python -c "import sys; print(sys.executable)"` shows a system Python instead of
the environment. Select the interpreter explicitly: `Ctrl+Shift+P` → "Python:
Select Interpreter" → **"Enter interpreter path…"** → the full path to
`envs/esp/python.exe`. Do not pick from the list if it still contains stale
entries. Then close *all* terminals and open a new one.

As a fallback that depends on no shell configuration at all, call the
environment's interpreter directly:

```powershell
& "C:\Users\<you>\miniforge3\envs\esp\python.exe" render_esp.py --prefix molecule
```

**`Unbekanntes Element 0.0000` / `Unknown element 0.0000` from `xyzToCube.py`.**
Fixed — the script now detects molfiles. If you see this on an older copy, the
structure file is an MDL molfile being parsed as xyz: molfiles list the
coordinates before the element symbol.

**Transparent surfaces show artefacts where atoms overlap.**
`set transparency_mode, 2` — without it PyMOL sorts transparent faces
incorrectly. The provided scripts set this already.

---

## 13. References

**Method / convention**

* F. A. Bulat, A. Toro-Labbé, T. Brinck, J. S. Murray, P. Politzer,
  *Quantitative analysis of molecular surfaces: areas, volumes, electrostatic
  potentials and average local ionization energies*,
  J. Mol. Model. **16** (2010) 1679–1691. [doi:10.1007/s00894-010-0692-x](https://doi.org/10.1007/s00894-010-0692-x)
* J. S. Murray, P. Politzer, *The electrostatic potential: an overview*,
  WIREs Comput. Mol. Sci. **1** (2011) 153–163. [doi:10.1002/wcms.19](https://doi.org/10.1002/wcms.19)

**Software**

* PyMOL — The PyMOL Molecular Graphics System, Schrödinger, LLC.
  Open-source build: <https://github.com/schrodinger/pymol-open-source>


---

*Practical Bioinformatics Project — Visualization of Molecular Electrostatic
Potentials.*
