================================================================================
tools/CreateTpTdFromSmiles.py - test data for the ESP workflow
================================================================================

WHAT THIS IS FOR
----------------

The workflow in scripts/ was developed on three halobenzenes: chloro-, bromo-
and iodobenzene. Structurally these are so similar that whole classes of bugs
cannot show up. Open questions were:

  * What happens when a molecule contains a group more negative than the
    halogen belt? In all three test molecules V_S,min sat on the halogen, so
    the "V_S,min" and "belt" lines always reported the same number. Whether
    the separation works at all was untested.

  * What happens with a molecule WITHOUT a halogen? There is no C-X axis for
    the orientation to align to, and the sigma-hole analysis has nothing to do.
    Does the script crash, or fall back cleanly?

  * Does xyzToCube.py also handle grids that do NOT come from Turbomole -
    different box dimensions, different origin, different point counts?

There is no public database of Turbomole pointval files. They are the
intermediate output of one particular calculation in one particular group;
nobody archives them. So we generate our own, aimed at exactly the cases we
are missing.


PLATFORM REQUIREMENT - READ FIRST
---------------------------------

PySCF is NOT available for Windows. conda-forge builds it for linux-64,
linux-aarch64, linux-ppc64le, macOS-64 and macOS-arm64 - there is no win-64
package, and there are no Windows wheels on PyPI either. On Windows,
'conda env create' will fail with:

    PackagesNotFoundError: The following packages are not available from
    current channels:
      - pyscf

This affects ONLY this tool. The actual workflow - xyzToCube.py,
render_esp.py, run_all.py - runs natively on Windows, macOS and Linux.

On Windows, use WSL (Windows Subsystem for Linux):

    1. In an Administrator PowerShell:
           wsl --install
       Reboot when asked, then set a user name and password.

    2. Inside the WSL shell, install Miniforge:
           curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
           bash Miniforge3-Linux-x86_64.sh
       Open a new shell afterwards.

    3. Your Windows files are mounted under /mnt/c, so the repository is at
           cd /mnt/c/Users/<you>/Desktop/.../Pymol_esp_visualization
       No copying needed - it is the same folder.

    4. Create the environment and run as described under USAGE below.

VS Code works with this directly; you do not have to drop to a bare console.
Install the "WSL" extension (ms-vscode-remote.remote-wsl), then either
F1 -> "WSL: Connect to WSL", or run 'code .' from a WSL shell inside the
project folder. The integrated terminal is then a Linux shell and
"Python: Select Interpreter" lists the conda environments living in WSL.

Note that writing across /mnt/c is slower than inside the Linux file system.
For the file sizes here (a few tens of MB) that is not a problem.


THE PIPELINE
------------

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

The result is a folder that looks exactly like a real dataset and can be
dropped straight into sandbox/.


WHY POINTVAL AND NOT CUBE DIRECTLY
----------------------------------

PySCF could write cube files directly. That is precisely what we do not want:
it would skip xyzToCube.py, and with it the Bohr/Angstrom conversion and the
reordering from x-fastest to z-fastest. That is the most error-prone step in
the whole chain. A test case that leaves it out tests the wrong thing.


PARAMETER CHOICES
-----------------

--spacing 0.25 Bohr (default)

    The real Turbomole data uses 0.12 Bohr; the cubes derived from it with
    --stride 2 use 0.24 Bohr. So 0.25 is essentially the same resolution -
    comparable, without blowing up the run time.

    Why not finer: evaluating the potential costs one integral over all basis
    functions per grid point. Halving the spacing multiplies the point count
    by eight. At 0.12 Bohr this would take hours instead of minutes, which is
    pointless for a functional test.

    Lower bound: render_esp.py determines the sigma-hole by casting rays and
    interpolating, not by picking grid points, so it is fairly tolerant of
    coarse grids. It still warns above 0.30 Bohr spacing, because the
    interpolated density itself smooths the isosurface and biases the value low
    by a few percent. 0.25 Bohr stays below that.

--margin 3.5 Angstrom (default)

    Distance from the outermost nucleus to the edge of the box. The
    rho = 0.001 isosurface sits roughly at the van der Waals radius, i.e. a
    good 2 Angstrom beyond the nuclei. 3.5 leaves headroom so the surface is
    not clipped - clipped isosurfaces show up in the images as straight edges.

    Bigger is not better: the margin applies in all three directions, so every
    extra Angstrom costs disproportionately more compute.

--basis def2-svp, --method hf (defaults)

    HF/def2-SVP is the classic level for ESP evaluation and it is fast.
    def2-SVP treats bromine all-electron; effective core potentials start at
    rubidium, so iodine uses one. B3LYP is available via --method but takes
    longer and changes nothing about the purpose of the test.

Block size of the potential evaluation

    No longer guessed, but derived from the memory footprint. int1e_grids
    returns an array of shape (points, nao, nao); with 193 basis functions
    that is 0.28 MB per single grid point. A fixed block of 20000 points tried
    to allocate 6 GB and was killed by the operating system. The block is now
    chosen so the intermediate stays below 400 MB.


SIGN OF THE POTENTIAL
---------------------

V(r) = sum_A Z_A/|R_A - r|  -  integral rho(r')/|r' - r| dr'

PySCF's int1e_grids returns the integrals with a POSITIVE sign, so the
electronic contribution has to be subtracted. Sanity check: far outside a
neutral molecule V must go to zero. For HCl at 40 Bohr, v_nuc = +0.4516 and
v_ele = +0.4512; the difference is +0.0004, while the sum would be +0.9028.
Get this backwards and the potential is positive everywhere - which is
immediately obvious, because no image has a single red region left.


LIMITATIONS - PLEASE READ
-------------------------

1. THE GEOMETRY COMES FROM A FORCE FIELD, not from a quantum chemical
   optimisation. MMFF94 gives usable structures, but bond lengths and angles
   deviate from an optimised geometry, which shifts the ESP values.

2. THE LEVEL OF THEORY IS NOT THE ONE USED FOR THE PROVIDED DATA. Method,
   basis set, and whether effective core potentials were used all affect
   V_S,min and V_S,max.

It follows that THESE NUMBERS DO NOT BELONG IN A TABLE NEXT TO VALUES FROM THE
PROVIDED TURBOMOLE DATA. They answer the question "does the pipeline run
correctly for this class of molecule", not "how large is the sigma-hole of
compound X". For the latter, every molecule would have to be optimised and
computed at the same level.

Incidentally this is the same footnote that is missing from the comparison of
the three halobenzenes: there, the maximum density value suggests bromine was
treated all-electron while iodine used an effective core potential. Worth
asking the supervisor.


BUILT-IN TEST CASES
-------------------

4-bromoacetophenone   CC(=O)c1ccc(Br)cc1
    Halogen AND carbonyl. The carbonyl oxygen is considerably more negative
    than the bromine belt, so V_S,min has to move onto the oxygen while the
    belt value is still measured at the bromine. This is the first case where
    the two lines show different numbers - for the halobenzenes they were
    always identical.

    Result: V_S,min = -41.0 kcal/mol on O3, belt = -10.5 kcal/mol at Br,
    sigma-hole = +15.2 kcal/mol. Note that the sigma-hole is almost twice that
    of bromobenzene (+7.9) and on par with iodobenzene (+15.5) - the acetyl
    group withdraws density from the ring and deepens the hole.

paracetamol           CC(=O)Nc1ccc(O)cc1
    No halogen. Tests whether the sigma-hole analysis is skipped cleanly and
    whether the orientation falls back to the principal axes instead of
    crashing.

    Result: V_S,min = -45.9 kcal/mol on O3, V_S,max = +56.1 kcal/mol on H18,
    no sigma-hole block, dash in the summary table, orientation sensible.


USAGE
-----

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
    python run_all.py --root ../sandbox

Expect roughly 13 seconds for the SCF and about 4 minutes for the grid
evaluation per molecule at 0.25 Bohr, with a progress display and an estimate
of the remaining time. Each of td.xyz and tp.xyz is around 60 MB.

The generated td.xyz/tp.xyz are large and excluded by .gitignore. That is
intentional: they are reproducible from the SMILES at any time with the
command above.


================================================================================
iso_sweep.py and stride_sweep.py - the parameter study
================================================================================

WHAT THESE ARE FOR
------------------

Two of the workflow's defaults are asserted in the background document rather
than obvious: rho = 0.001 for the isosurface, and stride 1 whenever a number
goes into a table. Both claims need measurements behind them, and those
measurements have to be repeatable - by the next reader, and by us after the
next change to render_esp.py. These two scripts are that measurement.

    iso_sweep.py       varies --iso at fixed resolution
                       -> table in section 4.1 of the background document

    stride_sweep.py    varies --stride at fixed isovalue
                       -> table in section 4.2 of the background document

They answer two different questions. The isovalue moves the physics: over
0.0005 .. 0.004 a.u. the sigma-hole of bromobenzene grows by a factor of 4.4,
so a sigma-hole value quoted without its isovalue means nothing. The stride
only costs accuracy: from 0.12 to 0.96 Bohr the same value falls by 17 %,
smoothly and always low, while the files shrink by a factor of 480. One
parameter has to be fixed by convention, the other can be traded against disk
space - and it is worth being able to show that, not just say it.

Both scripts import esp_statistics(), shell_points(), local_extrema() and
sigma_hole_interpolated() from ../scripts and call nothing of their own. That
is deliberate: a parameter study that reimplements the measurement proves
something about the study, not about the workflow. If render_esp.py changes,
these tables change with it.


WHAT THEY NEED
--------------

The normal 'esp' environment from environment.yml - numpy is enough. NOT the
esp-testdata environment above; there is no PySCF and no RDKit involved, so
these two run natively on Windows.

PyMOL is not required either, although render_esp.py is imported: that module
loads PyMOL only inside ensure_pymol(), which the measurement path never
reaches.

Input is td.cube and tp.cube, not the pointval files. Convert once with
xyzToCube.py, then both studies run on the result in seconds instead of parsing
1.25 GB of ASCII per pass.


USAGE
-----

    conda activate esp-pymol

    cd tools
    python iso_sweep.py    --folder ../sandbox/brombenzol
    python stride_sweep.py --folder ../sandbox/brombenzol

Each prints its table and writes a CSV next to the cube files
(iso_sweep_<folder>.csv, stride_sweep_<folder>.csv) carrying more columns than
the document shows - among them the point-based sigma-hole for comparison with
the ray-based one, which is the evidence for section 2.2.

Other sampling points, and another molecule:

    python iso_sweep.py    --isos 0.001 0.002 --folder ../sandbox/iodbenzol
    python stride_sweep.py --strides 1 2 4 --folder ../sandbox/chlorbenzol


ONE IMPLEMENTATION NOTE
-----------------------

stride_sweep.py does not reconvert the pointval files for every stride. It
decimates the full-resolution cube in memory with data[::N, ::N, ::N] and
scales the voxel vectors by N - which is exactly, line for line, what
write_cube() does when it is given --stride N. The results are therefore not
"comparable to" a real stride-N run; they are the same numbers.

The "cubes" column is computed from the cube format rather than measured,
because writing the stride-1 file only to read off its size would mean 200 MB
of disk traffic for one table cell. The data part of that formula is exact; the
header is estimated to within a few bytes, since the comment lines contain the
source file name. Checked against sandbox/brombenzol/td.cube: 210 865 313 bytes
predicted, 210 865 313 measured.
