# --------------------------------------------------------------
# esp.pml - ESP auf Elektronendichte-Isoflaeche
# Start:  pymol esp.pml       (oder in PyMOL:  @esp.pml)
# --------------------------------------------------------------

reinitialize

# 1) Struktur und Volumendaten laden
load brombenzol_aro_opti.mol, mol
load td.cube, dens
load tp.cube, esp

# 2) Molekuel als Staebchen
hide everything
show sticks, mol
set stick_radius, 0.3
color grey70, mol and elem C
util.cnc mol

# 3) Isoflaeche der Elektronendichte bei rho = 0.001 a.u.
#    (Politzer/Murray-Konvention fuer die "Molekueloberflaeche")
isosurface surf, dens, 0.001

# 4) Farbrampe fuer das ESP; Werte in Hartree/e (a.u.)
#    -0.035 .. 0.035 a.u.  entspricht -92 .. 92 kJ/(mol*e)
ramp_new espramp, esp, [-0.035, 0, 0.035], [red, white, blue]

# 5) ESP auf die Oberflaeche mappen
set surface_color, espramp, surf
set surface_quality, 1

#    Transparenz: 0 = opak (kraftigste Farben, Staebchen unsichtbar),
#    0.15 = Standard (Molekuelgeruest scheint durch),
#    ab ca. 0.3 wird es unleserlich, weil man durch das ganze Molekuel schaut.
set transparency, 0.15
set transparency_mode, 2
set two_sided_lighting, on

# 6) Darstellung / Rendering
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 0
set specular, 0.2
set ambient, 0.15
orient mol
zoom mol, 2.0

# 7) Hochaufloesendes Bild
# ray 2400, 1800
# png esp.png, dpi=300
