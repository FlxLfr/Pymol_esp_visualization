# `tools/` creating additional test data, creating the reference dataset and conducting the parameter study

Four scripts that sit **outside** the rendering workflow. None of them is called
by `run_all.py`, and nothing in `scripts/` imports them.

| Script | Purpose | Documented in |
|---|---|---|
| `CreateTpTdFromSmiles.py` | generates Turbomole-format test data from SMILES | §1 of this file |
| `make_reference.py` | builds the committed self-test dataset | §2 of this file |
| `iso_sweep.py` | varies the isovalue → ProjectElaboration.pdf §4.1 | §3 of this file |
| `stride_sweep.py` | varies the grid stride → ProjectElaboration.pdf §4.2 | §3 of this file |

---

## 1  `CreateTpTdFromSmiles.py`, creating edge cases data for the ESP workflow

### What this is for

The workflow in `scripts/` was developed on three halobenzenes and some pyridine
derivatives. Structurally these are so similar that whole classes of bugs cannot
show up.

There is no public database of Turbomole `pointval` files. They are the
intermediate output of one particular calculation in one particular group;
nobody archives them. So the need to generate our own with the help of this script arised.

### Platform requirement!!

PySCF is **not** available for Windows. conda-forge builds it for linux-64,
linux-aarch64, linux-ppc64le, macOS-64 and macOS-arm64. There is no win-64
package, and there are no Windows wheels on PyPI either. On Windows,
`conda env create` will fail with:

```
PackagesNotFoundError: The following packages are not available from
current channels:
  - pyscf
```

This affects **only this tool**. The actual esp visualization workflow runs natively on Windows, macOS and Linux.

On Windows, use WSL (Windows Subsystem for Linux):

1. In an Administrator PowerShell:

   ```powershell
   wsl --install
   ```

   Reboot when asked, then set a user name and password.

2. Inside the WSL shell, install Miniforge:

   ```bash
   curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh
   ```

   Open a new shell afterwards.

3. Your Windows files are mounted under `/mnt/c`, so the repository is at

   ```bash
   cd /mnt/c/Users/<you>/Desktop/.../Pymol_esp_visualization
   ```

4. Create the environment and run as described under [Usage](#usage) below.

Note that writing across `/mnt/c` is slower than inside the Linux file system.
For the file sizes here (a few tens of MB) that is not a problem.

### The pipeline

```
SMILES
  |  RDKit: ETKDGv3 embedding + MMFF94 optimisation
  v
3D geometry  ---> <name>.mol
  |  PySCF: SCF (default HF/def2-SVP)
  v
density matrix
  |  evaluated on a regular grid
  v
rho(r) ---> td.xyz        V(r) ---> tp.xyz     (Turbomole pointval format)
```

The result is a folder that looks exactly like a real dataset and can be dropped
straight into `sandbox/`.


### Options

| Option | Default | Effect |
|---|---|---|
| `--smiles` | - | SMILES of a single molecule; use together with `--name` |
| `--name` | - | folder and file name for that molecule |
| `--preset` | off | build the four built-in test cases instead: 4-bromoacetophenone, paracetamol, haloperidol and triazolam |
| `--outdir` | `../sandbox` | root directory the molecule folders are written into |
| `--spacing` | `0.25` | grid spacing in Bohr |
| `--margin` | `3.5` | margin around the molecule, in Ångström |
| `--basis` | `def2-svp` | basis set handed to PySCF |
| `--method` | `hf` | SCF method, `hf` or `b3lyp` |

Either `--preset` or `--smiles` together with `--name` has to be given; without
one of the two the script stops and says so. Why the defaults are what they are
is the next section.

### Parameter choices

**`--spacing 0.25` Bohr (default).** The real Turbomole data uses 0.12 Bohr; the
cubes derived from it with `--stride 2` use 0.24 Bohr. So 0.25 is essentially the
same resolution, comparable, good enough for testing, without blowing up the run time.

Why not finer: evaluating the potential costs one integral over all basis
functions per grid point. Halving the spacing multiplies the point count by
eight. At 0.12 Bohr this would take hours instead of minutes, which is pointless
for a functional test.

**`--margin 3.5` Ångström (default).** Distance from the outermost nucleus to the
edge of the box. The ρ = 0.001 isosurface sits roughly at the van der Waals
radius, i.e. a good 2 Å beyond the nuclei. 3.5 leaves headroom so the surface is
not clipped — clipped isosurfaces show up in the images as straight edges.

Bigger is not better: the margin applies in all three directions, so every extra
Ångström costs disproportionately more computing time.

**`--basis def2-svp`, `--method hf` (defaults).** HF/def2-SVP is the classic level
for ESP evaluation and it is fast. def2-SVP treats bromine all-electron;
effective core potentials start at rubidium, so iodine uses one. B3LYP is
available via `--method` but takes longer and changes nothing about the purpose
of the test.


### Limitations

1. **The geometry comes from a force field**, not from a quantum chemical
   optimisation. MMFF94 gives usable structures, but bond lengths and angles
   deviate from an optimised geometry, which shifts the ESP values.
2. **The level of theory is not the one used for the provided data.** Method,
   basis set and whether effective core potentials were used all affect V_S,min
   and V_S,max.

> **These numbers do not belong in a table next to values from the provided
> Turbomole data.** They answer the question "does the pipeline run correctly for
> this other class of molecule", not "how large is the σ-hole of compound X". For the
> latter, every molecule would have to be optimised and computed at the same
> level.

### Built-in test cases

**4-bromoacetophenone**: `CC(=O)c1ccc(Br)cc1`

Halogen **and** carbonyl. The carbonyl oxygen is considerably more negative than
the bromine belt, so V_S,min has to move onto the oxygen while the belt value is
still measured at the bromine. This is the first case where the two lines show
different numbers: for the halobenzenes they were always identical.

Result: V_S,min = −41.0 kcal/mol on O3, belt = −10.5 kcal/mol at Br, σ-hole =
+15.2 kcal/mol. That is about 1.5 times the +10.2 of bromobenzene and close to
the +16.1 of iodobenzene. The acetyl group withdraws density from the ring and
deepens the hole.

**Paracetamol**: `CC(=O)Nc1ccc(O)cc1`

No halogen. Tests whether the σ-hole analysis is skipped cleanly and whether the
orientation falls back to the principal axes instead of crashing.

Result: V_S,min = −45.9 kcal/mol on O3, V_S,max = +56.1 kcal/mol on H18, no
σ-hole block, dash in the summary table, orientation sensible.

**Triazolam (Halcion)**: `CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4Cl`

2 halogens and super complex 3D-molecule structure, here two errors occurred with the implementation, 
the first halogen was taken for the axial view picture and not the one with the biggest sigma hole and the
rays to calculate the sigma hole, left the halogen and measured the surface of the neighbouring methyl group (see ProjectElaboration.pdf §5.2).
Both defects have been subsequently fixed and Halcion was added to the results, as it was crucial for the development.

Result: V_S,min = -52.89 kcal/mol on N4, V_S,max = +34.09 kcal/mol on H27
σ-hole (Cl21) = +0.01709 a.u. (+10.73 kcal/(mol*e))  [interpolated]
belt (Cl21) = -0.01836 a.u. (-11.52 kcal/(mol*e))
σ-hole (Cl11) = +0.01440 a.u. (+9.03 kcal/(mol*e))  [interpolated]
belt (Cl11) = -0.01095 a.u. (-6.87 kcal/(mol*e))


**Haloperidol**: `O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1`

Used to confirm the fixes, implemented by the findings of Triazolam, also 2 halogens are available and the sigma view showed the halogen with the bigger potential value 

Result: V_S,min = -42.42 kcal/mol  auf O19, V_S,max = +43.98 kcal/mol on H39
σ-hole (Cl13) = +0.00551 a.u. (+3.46 kcal/(mol*e))  [interpolated]
belt (Cl13) = -0.01925 a.u. (-12.08 kcal/(mol*e))
σ-hole (F26) = -0.02631 a.u. (-16.51 kcal/(mol*e))  [interpolated]
belt (F26) = -0.02142 a.u. (-13.44 kcal/(mol*e))

### Usage

```bash
# Linux / macOS, or inside WSL on Windows
conda env create -f tools/environment-testdata.yml
conda activate esp-testdata

cd tools
python CreateTpTdFromSmiles.py --preset

# or a molecule of your own
python CreateTpTdFromSmiles.py --smiles "O=C(N)c1ccccc1" --name benzamide

# quick check that the chain works, four seconds instead of four minutes
python CreateTpTdFromSmiles.py --smiles "ClC" --name chloromethane \
                               --spacing 0.4 --margin 3.0

# then continue as usual
cd ../scripts
python run_all.py --root ../sandbox --only chloromethane
```

---

## 2  `make_reference.py`

### What this is for

`python run_all.py` without arguments runs on `reference/brombenzol/`. That is
the smoke test after a fresh clone: it converts, renders, and writes to separate
`_check` names (`images_check/`, `summary_check_<time>_<date>.csv`), so the
committed reference output stays untouched.

The test answers exactly one question: **does the installation run, and do the
documented numbers come out?** It does not answer "is this an accurate V_S value"

The expected result on this dataset:

```
V_S,min = -0.01863   V_S,max = +0.03070 a.u.  -> colour scale +/- 0.0350
sigma hole (Br12) = +0.01528 a.u. (+9.59 kcal/(mol*e))
```

If those come out, unit conversion, index reordering, the shell band, the ray
search and the colour-scale rule all work. If they do not, the problem is the
installation or a change to the scripts.

### How the dataset is built

```bash
python tools/make_reference.py sandbox/brombenzol --name brombenzol
```

It writes a decimated `td.xyz` / `tp.xyz` pair in Turbomole `pointval` format to
`reference/<name>/` and copies the structure file along — without it `run_all.py`
will not recognise the folder.

| Option | Default | Effect |
|---|---|---|
| `source` | *required* | molecule folder with `td`/`tp` and a structure file |
| `--name` | the folder name | name under `reference/` |
| `--outdir` | the repository's `reference/` | write somewhere else |
| `--stride N` | `5` | keep every N-th grid point per axis |
| `--margin` | `2.5` | margin around the isosurface, in Bohr |
| `--iso` | `0.001` | the isovalue that must stay fully contained |

### Two steps, both necessary

**Cropping.** The full grid is a 30 Bohr box and the molecule with its
ρ = 0.001 surface fills only the middle of it. The crop is derived from the
**density**: the bounding box of all points with ρ > iso/2, plus `--margin`. That
guarantees the isosurface is fully contained and the images are not cut off at
the edge.

**Decimating.** Then only every `--stride`-th point per axis.

Together the two bring 1.25 GB down to 1.7 MB per `pointval` file.

### Why `pointval` and not cube

The smoke test is meant to exercise the whole chain, including `xyzToCube.py`,
unit conversion and index reordering are the two steps most likely to break.
Shipping ready-made cubes would skip exactly those.

### What it costs

0.60 Bohr is five times coarser than the production grids at 0.12 Bohr, and that
is measurable:

| | reference, 0.60 Bohr | production, 0.12 Bohr | deviation |
|---|---|---|---|
| V_S,min | −0.01863 | −0.01882 | 1.0 % |
| V_S,max | +0.03070 | +0.03154 | 2.7 % |
| σ-hole | +0.01528 | +0.01629 | **6.2 %** |

The V_S values hold to within one to three per cent; the σ-hole falls short by
six. That is the smooth, one-sided degradation the stride study in section 4.2 of
ProjectElaboration.pdf describes: the interpolated density smooths the isosurface
outwards and the potential is less positive there. Broad flat features barely
notice, sharp ones suffer.

The reason for the coarse grid is git, not physics. The dataset has to be
committed so that a fresh clone can test itself, and 1.7 MB per file is
committable where 200 MB is not.

> **Do not quote a V_S or σ-hole value from the reference dataset.** It exists to
> prove the chain works. Numbers that go into a table come from the full-resolution
> grids under `sandbox/` or `results/`.

---

## 3  `iso_sweep.py` and `stride_sweep.py`

### What these are for

Two of the workflow's defaults are asserted in ProjectElaboration.pdf rather than
obvious: ρ = 0.001 for the isosurface and stride 1 whenever a number goes into a
table. Both claims need measurements behind them and those measurements have to
be repeatable. These two scripts are that measurement.

| Script | Varies | Feeds |
|---|---|---|
| `iso_sweep.py` | `--iso` at fixed resolution | the table in section 4.1 of ProjectElaboration.pdf |
| `stride_sweep.py` | `--stride` at fixed isovalue | the table in section 4.2 of ProjectElaboration.pdf |

They answer two different questions. The isovalue moves the physics: over
0.0005 … 0.004 a.u. the σ-hole of bromobenzene grows by a factor of 4.4, so a
σ-hole value quoted without its isovalue means nothing. The stride only costs
accuracy: from 0.12 to 0.96 Bohr the same value falls by 17 %, smoothly and
always low, while the files shrink by a factor of 480. One parameter has to be
fixed by convention, the other can be traded against disk space.

Both scripts import `esp_statistics()`, `shell_points()`, `local_extrema()` and
`sigma_hole_interpolated()` from `../scripts` and call nothing of their own. That
is deliberate: a parameter study that reimplements the measurement proves
something about the study, not about the workflow. If `render_esp.py` changes,
these tables change with it.

### What they need

The normal `esp` environment from `environment.yml` is enough. **Not** the
`esp-testdata` environment above; there is no PySCF and no RDKit involved, so
these two run natively on Windows.

PyMOL is not required either, although `render_esp.py` is imported: that module
loads PyMOL only inside `ensure_pymol()`, which the measurement path never
reaches.

Input is `td.cube` and `tp.cube`, not the `pointval` files. Convert once with
`xyzToCube.py`, then both studies run on the result in seconds instead of parsing
1.25 GB of ASCII per pass.

### Usage

```bash
conda activate esp-pymol

cd tools
python iso_sweep.py    --folder ../sandbox/brombenzol
python stride_sweep.py --folder ../sandbox/brombenzol
```

Each prints its table and writes a CSV next to the cube files
(`iso_sweep_<folder>.csv`, `stride_sweep_<folder>.csv`) carrying more columns than
the document shows. Among them the point-based σ-hole for comparison with the
ray-based one, which is the evidence for section 4.3 of ProjectElaboration.pdf.

Other sampling points and another molecule:

```bash
python iso_sweep.py    --isos 0.001 0.002 --folder ../sandbox/iodbenzol
python stride_sweep.py --strides 1 2 4 --folder ../sandbox/chlorbenzol
```

### One implementation note

`stride_sweep.py` does not reconvert the `pointval` files for every stride. It
decimates the full-resolution cube in memory with `data[::N, ::N, ::N]` and scales
the voxel vectors by N, which is exactly, line for line, what `write_cube()` does
when it is given `--stride N`. The results are therefore not "comparable to" a
real stride-N run; they are the same numbers.

