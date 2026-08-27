# --------------------------------------------------------------
# esp.pml - ESP on the electron density isosurface
# Start:  pymol esp.pml       (or inside PyMOL:  @esp.pml)
# --------------------------------------------------------------

reinitialize

# 1) load the structure and the volumetric data
load brombenzol_aro_opti.mol, mol
load td.cube, dens
load tp.cube, esp

# 2) the molecule as sticks
hide everything
show sticks, mol
set stick_radius, 0.3
color grey70, mol and elem C
util.cnc mol

# 3) isosurface of the electron density at rho = 0.001 a.u.
#    (the Politzer/Murray convention for the "molecular surface")
isosurface surf, dens, 0.001

# 4) colour ramp for the ESP; values in Hartree/e (a.u.)
#    -0.035 .. 0.035 a.u.  equals -92 .. 92 kJ/(mol*e)
ramp_new espramp, esp, [-0.035, 0, 0.035], [red, white, blue]

# 5) map the ESP onto the surface
set surface_color, espramp, surf
set surface_quality, 1

#    transparency: 0 = opaque (strongest colours, sticks invisible),
#    0.15 = the default (the skeleton shows through),
#    from about 0.3 on it becomes unreadable, because you look through the
#    whole molecule.
set transparency, 0.15
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
