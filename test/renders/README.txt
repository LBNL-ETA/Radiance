# RCSid $Id: README.txt,v 1.2 2018/12/04 00:32:20 greg Exp $

ref/	Subdirectory containing reference outputs for comparison

Makefile	Scene testing dependencies

README.txt	This file
test.txt	Some text for use in *text primitives

render.opt	A few global rendering options used by most rad input files

fish.vf	View parameters for standard fisheye view away from window
inside.vf	View parameters for default interior view towards window

combined.rif	A combined scene file including every primitive
dielectric.rif	Test of dielectric material in window and mesh instance
glass.rif	Test of glass pane in window and material mixture
inst.rif	Test of octree instances and ashik2 material
mesh.rif	Test of Radiance triangle mesh instances and view types
mirror.rif	Test of mirror (virtual sources) and antimatter
mist.rif	Test of mist (participating medium) in combination with spotlight
mixtex.rif	Test of BSDF types and photon-mapping
patterns.rif	Test of various pattern types and alias behavior
prism1.rif	Test of prism1 type (single virtual source)
prism2.rif	Test of prism2 type (dual virtual sources)
tfunc.rif	Test of BRTDfunc, transfunc and mkillum
trans.rif	Test of anisotropic types, plasfunc and metfunc, mixtures
trans2.rif	Test of trans2 and mkillum (again)

basic.mat	Basic materials used by all models
chrome.mat	Chrome material used in mesh model
gold.mat	Gold material used in prism1 and prism2 models
mixtex.mat	Nine test materials used in mixtex model
patterns.mat	Nine test materials and aliases used in patterns model

antimatter_portal.rad	Cut-away hole seen in mirror model
ball_in_cube.rad	Example of dielectric and interface materials
ballcompare.rad	Combo test of *func materials, mixtures, and text for trans model
blinds.rad	BSDF from XML file with geometry placed by pkgBSDF
bubble_in_cube.rad	Example of colored dielectric and interface
closed_end.rad	Cap for far end of space when there is no window
combined_scene.rad	Combination of all our test scenes without front caps
constellation.rad	A set of ten simple light bulbs ina  circle
dielectric_pane.rad	Dual-surface dielectric window for testing
diorama_walls.rad	The walls of our deep test room, minus both ends
disks.rad	A generic set of nine disks for material testing
front_cap.rad	Wall away from window (used except in combined scene)
glass_ill.rad	Recorded output from mkillum run of glass window
glass_pane.rad	Bluish glass pane for window
glowbulb.rad	Test of glow material with non-zero influence radius
gymbal.rad	Example of BRTDfunc type, cylinders and tubes
illum_glass.rad	Input used to create glass_ill.rad
illum_tfunc.rad	Mkillum input used in tfunc model
illum_trans2.rad	Mkillum input used in trans2 model
mirror.rad	Example of mirror on wall
mist.rad	Example of mist (participating medium) inside test space
porsches.rad	Two instances of porsche.octf with original material and ashik2
prism1.rad	Example of rectangular prism1 pane in window
prism2.rad	Example of rectangular prism2 pane in window
rect_fixture.rad	A single rectangular light fixture with distribution
rect_opening.rad	Rectangular opening in window wall
saucer.rad	A shaped disk for material testing
sawtooth.rad	Object with sawtooth profile for genBSDF input (not referenced)
spotcones.rad	Cones enclosing spotlight regions for testing mist type
spotlights.rad	Color spotlights
sunset_sky.rad	Captured sky with disk covering sun position
torus.rad	Psychedelic donut testing material mixtures and glow
trans_pane.rad	A yellowish diffusing test of trans type
vase.rad	Instance of vase triangle mesh
woman.rad	Instance of "woman" triangle mesh

vase.rtm	Vase triangle mesh with local texture map
woman.rtm	Woman triangle mesh with local texture map
porsche.octf	Frozen octree of Porsche with materials

aniso.cal	Ward-Geisler-Moroder-Duer anisotropic reflection model
bumpypat.cal	Used by texdata type in patterns model
climit.cal	Standard color handling for local texture map used by vase.rad
diskcoords.cal	Coordinates used for Shirley-Chiu BSDF data
fisheye.cal	Mapping for captured sunset image
maxang.cal	Used for transdata material in mixtex.mat
prism.cal	Standard angle calculation for prismatic glazing (prism1 & prism2)

flower.hdr	Image of flower used in various *pict types
sunset.hdr	Captured sunset for environment map out window
vase.hdr	Vase local texture map

flower.dat	Flower.hdr converted to a 10x10 gray data file
glass_illB.dat	Distribution computed by mkillum for glass_ill.rad
glass_illG.dat	Distribution computed by mkillum for glass_ill.rad
glass_illR.dat	Distribution computed by mkillum for glass_ill.rad
rect_fixture.dat	Rectangular light fixture output distribution
sawtooth.dat	BRDF of sawtooth profile material converted to a data file
tcutoff.dat	Used for transdata material in mixtex.mat
