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
2. [Software: what was evaluated and what was chosen](#2-software-what-was-evaluated-and-what-was-chosen)
3. [Installation](#3-installation)
4. [Input files and formats](#4-input-files-and-formats)
5. [Step 1 — Convert the grids to cube](#5-step-1--convert-the-grids-to-cube)
6. [Step 2 — Look at it interactively](#6-step-2--look-at-it-interactively)
7. [Step 3 — Render the standard image set](#7-step-3--render-the-standard-image-set)
8. [Step 4 — Several molecules at once](#8-step-4--several-molecules-at-once)
9. [Choosing the colour scale](#9-choosing-the-colour-scale)
10. [Repository layout](#10-repository-layout)
11. [Troubleshooting](#11-troubleshooting)
12. [References](#12-references)

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
[§7](#which-number-describes-the-σ-hole).

---

## 2. Software: what was evaluated and what was chosen

| Program | Reads cube | Surface + ESP mapping | Colour scale control | Scriptable | Licence |
|---|---|---|---|---|---|
| **PyMOL (open source)** | yes | `isosurface` + `ramp_new` | full, arbitrary ramps | **Python API + `.pml`** | free, open source |
| VMD | yes | Isosurface + Colorvolume | full | Tcl / Python | free, academic |
| UCSF ChimeraX | yes | `surface` + `color electrostatic` | full | command scripts + Python | free, non-commercial |
| Avogadro 2 | yes | yes | limited | limited | free, open source |
| GaussView | yes | yes | limited | **no** | commercial |
| Chemcraft | yes | yes | moderate | **no** | commercial |
| Multiwfn | (generates them) | quantitative analysis | n/a | shell scripting | free |

**PyMOL was chosen** because it is the only candidate that combines all four
requirements of this project: it reads Gaussian cube directly, it maps one volume
onto an isosurface of another, its entire state is reachable from Python — so the
figure settings live in a file under version control rather than in a sequence of
mouse clicks — and it is free and open source, so the SOP can specify the exact
installation as one command instead of a licence request.

GaussView and Chemcraft were ruled out on the automation criterion alone: neither
can be driven from a script, so neither can guarantee that two molecules were
rendered with identical settings.

**Multiwfn** is worth adding to the toolchain later if quantitative surface
descriptors are needed (V<sub>S,min</sub>/V<sub>S,max</sub> statistics, surface
areas, σ-hole magnitudes). This workflow already reports V<sub>S,min</sub> and
V<sub>S,max</sub>, but Multiwfn goes considerably further.

---

## 3. Installation

### 3.1 Prerequisite: a working conda

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
> [Troubleshooting](#11-troubleshooting) for how to get out of that.

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

### 3.2 Create the environment

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

### 3.3 Verify

```bash
python -c "import pymol, numpy, matplotlib; print('ok')"
python -c "import sys; print(sys.executable)"
```

The second line must point *inside* the `esp` environment. If it points at a
system Python instead, the environment was not activated — see
[Troubleshooting](#11-troubleshooting).

---

## 4. Input files and formats

Per molecule, in one folder:

| File | Content | Units |
|---|---|---|
| `td.xyz` | Turbomole `pointval` **total density** grid | Bohr |
| `tp.xyz` | Turbomole `pointval` **total potential** (ESP) grid | Bohr |
| `*.mol` / `*.sdf` / `*.pdb` / `*.xyz` | molecular structure | Å (default) |

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

Both `xyzToCube.py` and `render_esp.py` accept the same four formats:

| Format | Notes |
|---|---|
| `.xyz` | `Symbol x y z`. With or without the leading atom-count and comment lines — a bare coordinate list is accepted. Unit set by `--struct-unit` (default Å). |
| `.mol` | MDL molfile, V2000 and V3000. **Coordinates come *before* the element symbol**, the reverse of xyz. Always Å. |
| `.sdf` | SD-file; only the first record (up to `$$$$`) is read. Always Å. |
| `.pdb` | `ATOM`/`HETATM` records. Element from columns 77–78, otherwise derived from the atom name. Always Å. |

`--struct-unit` applies to `.xyz` only. Molfile and PDB are in Ångström by
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

## 5. Step 1 — Convert the grids to cube

```bash
cd scripts
python xyzToCube.py --struct ../path/to/molecule.mol ../path/to/td.xyz ../path/to/tp.xyz --pymol
```

`--struct` takes `.xyz`, `.mol`, `.sdf` or `.pdb` — see
[§4](#accepted-structure-formats).

This writes `td.cube`, `tp.cube` and (with `--pymol`) a ready-to-use `esp.pml`
next to the input.

Useful options:

| Option | Effect |
|---|---|
| `--stride 2` | keep every 2nd grid point in each direction → **8× smaller** files |
| `--struct-unit bohr` | structure file is already in Bohr |
| `--esp-range 0.035` | colour range written into the generated `esp.pml` |
| `--transparency 0` | opaque surface in the generated `esp.pml` |
| `--outdir DIR` | write the cubes somewhere else |

**On `--stride`.** For bromobenzene, decimating the 251³ grid to 126³ leaves
V<sub>S,min</sub> unchanged and shifts V<sub>S,max</sub> by 0.9 % (+19.58 vs.
+19.80 kcal/(mol·e)). The rendered images are indistinguishable, the files shrink
from 201 MB to 26 MB and PyMOL becomes noticeably faster.

**The σ-hole is more sensitive.** It is evaluated by ray casting with
interpolation rather than read off grid points, which makes it far more robust,
but a coarse grid still smooths the isosurface and biases the value low by a few
percent (see [§7](#which-number-describes-the-σ-hole)). So: `--stride 2` while
you are exploring and for the images, full resolution whenever a σ-hole value
goes into a table.

What the converter takes care of, which is where hand-rolled conversions usually
go wrong:

* **Index order.** Turbomole varies *x* fastest, the cube format varies *z*
  fastest. Without reordering you get a transposed, mirrored molecule.
* **Units.** The grid is in Bohr, the structure file is normally in Å (factor
  1.8897). Get this wrong and the molecule floats outside its own surface.

---

## 6. Step 2 — Look at it interactively

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

## 7. Step 3 — Render the standard image set

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
no manual rotation is needed or wanted.

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

`--stride 2` is fine for the images and for V<sub>S,min</sub> /
V<sub>S,max</sub> — those change by about 1 %. For σ-hole values that go into a
table, use the full grid.

Options worth knowing:

| Option | Default | Effect |
|---|---|---|
| `--esp-range` | `auto` | `auto` or a fixed value in a.u. — see [§9](#9-choosing-the-colour-scale) |
| `--transparency` | `0.15` | `0` = opaque, strongest colours; above ~0.3 the profile views become unreadable |
| `--backgrounds white black` | `white` | render each view on both backgrounds |
| `--views pi sigma` | all three | subset of views |
| `--width / --height / --dpi` | 2000 / 1600 / 300 | image size |
| `--iso` | `0.001` | density isovalue |
| `--buffer` | `2.4` | margin around the molecule, Å |

**Do not use the PyMOL launcher unless you have to.** `python render_esp.py …`
loads no `pymolrc`, so nobody's personal start-up file can silently change a
setting. If you do use the launcher, both the `--` separator and `-k` are
required:

```bash
pymol -ckq render_esp.py -- --prefix molecule
```

---

## 8. Step 4 — Several molecules at once

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

## 9. Choosing the colour scale

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

---

## 10. Repository layout

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
├── docs/                         exported PDF of this SOP
└── sandbox/                      your own data and experiments, not tracked
```

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

## 11. Troubleshooting

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

## 12. References

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
* W. Humphrey, A. Dalke, K. Schulten, *VMD — Visual Molecular Dynamics*,
  J. Mol. Graph. **14** (1996) 33–38. [doi:10.1016/0263-7855(96)00018-5](https://doi.org/10.1016/0263-7855(96)00018-5)
* E. C. Meng et al., *UCSF ChimeraX: Tools for structure building and analysis*,
  Protein Sci. **32** (2023) e4792. [doi:10.1002/pro.4792](https://doi.org/10.1002/pro.4792)
* T. Lu, F. Chen, *Multiwfn: A multifunctional wavefunction analyzer*,
  J. Comput. Chem. **33** (2012) 580–592. [doi:10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)
* T. Lu, *A comprehensive electron wavefunction analysis toolbox for chemists,
  Multiwfn*, J. Chem. Phys. **161** (2024) 082503. [doi:10.1063/5.0216272](https://doi.org/10.1063/5.0216272)
* M. D. Hanwell et al., *Avogadro: an advanced semantic chemical editor,
  visualization, and analysis platform*, J. Cheminform. **4** (2012) 17.
  [doi:10.1186/1758-2946-4-17](https://doi.org/10.1186/1758-2946-4-17)
* G. Schaftenaar, J. H. Noordik, *Molden: a pre- and post-processing program for
  molecular and electronic structures*, J. Comput.-Aided Mol. Des. **14** (2000)
  123–134. [doi:10.1023/A:1008193805436](https://doi.org/10.1023/A:1008193805436)

---

*Practical Bioinformatics Project — Visualization of Molecular Electrostatic
Potentials.*
