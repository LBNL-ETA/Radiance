# RCSid $Id: Makefile,v 1.5 2018/11/20 01:49:16 greg Exp $
#
# Render and primitive type testing
#

# Image reduction for comparisons
RDU_PFILT = pfilt -1 -r 1 -x 128 -y 128 -pa 1

# Image comparison command
IMG_CMP = radcompare -rms 0.08 -max 1.5

# Default target is to test everything
all:	test-dielectric-def test-dielectric-fish test-glass-def \
test-glass-fish test-glass-up test-inst-def test-inst-fish \
test-mesh-def test-mesh-cyl test-mirror-fish test-mist-def \
test-tfunc-def test-tfunc-fish test-trans-def test-trans-fish \
test-trans2-def test-trans2-fish

clean:
	rm -f *.oct *.amb *_ill.dat *_*.hdr *.unf combined.rad

# Special test of xform

test-xform:	combined_scene.rad
	xform -f combined_scene.rad | grep -v '^[ 	]*#' > combined.rad
	radcompare ref/combined.rad combined.rad

### Octree for dielectric ###

dielectric.oct:	dielectric.rif
	rad -v 0 dielectric.rif

### Reference and test for dielectric view   def ###

test-dielectric-def:	ref/dielectric_def.hdr dielectric_def.hdr
	$(RDU_PFILT) dielectric_def.hdr | $(IMG_CMP) ref/dielectric_def.hdr -

ref/dielectric_def.hdr:
	make dielectric_def.hdr
	$(RDU_PFILT) dielectric_def.hdr > ref/dielectric_def.hdr

dielectric_def.hdr:	dielectric.oct
	rad -v   def dielectric.rif

### End dielectric-def   tests

### Reference and test for dielectric view  fish ###

test-dielectric-fish:	ref/dielectric_fish.hdr dielectric_fish.hdr
	$(RDU_PFILT) dielectric_fish.hdr | $(IMG_CMP) ref/dielectric_fish.hdr -

ref/dielectric_fish.hdr:
	make dielectric_fish.hdr
	$(RDU_PFILT) dielectric_fish.hdr > ref/dielectric_fish.hdr

dielectric_fish.hdr:	dielectric.oct
	rad -v  fish dielectric.rif

### End dielectric-fish  tests

### Octree for   glass ###

glass.oct:	glass.rif
	rad -v 0   glass.rif

### Reference and test for   glass view   def ###

test-glass-def:	ref/glass_def.hdr   glass_def.hdr
	$(RDU_PFILT)   glass_def.hdr | $(IMG_CMP) ref/glass_def.hdr -

ref/glass_def.hdr:
	make   glass_def.hdr
	$(RDU_PFILT)   glass_def.hdr > ref/glass_def.hdr

glass_def.hdr:	glass.oct
	rad -v   def   glass.rif

### End   glass-def   tests

### Reference and test for   glass view  fish ###

test-glass-fish:	ref/glass_fish.hdr   glass_fish.hdr
	$(RDU_PFILT)   glass_fish.hdr | $(IMG_CMP) ref/glass_fish.hdr -

ref/glass_fish.hdr:
	make   glass_fish.hdr
	$(RDU_PFILT)   glass_fish.hdr > ref/glass_fish.hdr

glass_fish.hdr:	glass.oct
	rad -v  fish   glass.rif

### End   glass-fish  tests

### Reference and test for   glass view    up ###

test-glass-up:	ref/glass_up.hdr   glass_up.hdr
	$(RDU_PFILT)   glass_up.hdr | $(IMG_CMP) ref/glass_up.hdr -

ref/glass_up.hdr:
	make   glass_up.hdr
	$(RDU_PFILT)   glass_up.hdr > ref/glass_up.hdr

glass_up.hdr:	glass.oct
	rad -v    up   glass.rif

### End   glass-up    tests

### Octree for    inst ###

inst.oct:	inst.rif
	rad -v 0    inst.rif

### Reference and test for    inst view   def ###

test-inst-def:	ref/inst_def.hdr    inst_def.hdr
	$(RDU_PFILT)    inst_def.hdr | $(IMG_CMP) ref/inst_def.hdr -

ref/inst_def.hdr:
	make    inst_def.hdr
	$(RDU_PFILT)    inst_def.hdr > ref/inst_def.hdr

inst_def.hdr:	inst.oct
	rad -v   def    inst.rif

### End    inst-def   tests

### Reference and test for    inst view  fish ###

test-inst-fish:	ref/inst_fish.hdr    inst_fish.hdr
	$(RDU_PFILT)    inst_fish.hdr | $(IMG_CMP) ref/inst_fish.hdr -

ref/inst_fish.hdr:
	make    inst_fish.hdr
	$(RDU_PFILT)    inst_fish.hdr > ref/inst_fish.hdr

inst_fish.hdr:	inst.oct
	rad -v  fish    inst.rif

### End    inst-fish  tests

### Octree for    mesh ###

mesh.oct:	mesh.rif
	rad -v 0    mesh.rif

### Reference and test for    mesh view   def ###

test-mesh-def:	ref/mesh_def.hdr    mesh_def.hdr
	$(RDU_PFILT)    mesh_def.hdr | $(IMG_CMP) ref/mesh_def.hdr -

ref/mesh_def.hdr:
	make    mesh_def.hdr
	$(RDU_PFILT)    mesh_def.hdr > ref/mesh_def.hdr

mesh_def.hdr:	mesh.oct
	rad -v   def    mesh.rif

### End    mesh-def   tests

### Reference and test for    mesh view   cyl ###

test-mesh-cyl:	ref/mesh_cyl.hdr    mesh_cyl.hdr
	$(RDU_PFILT)    mesh_cyl.hdr | $(IMG_CMP) ref/mesh_cyl.hdr -

ref/mesh_cyl.hdr:
	make    mesh_cyl.hdr
	$(RDU_PFILT)    mesh_cyl.hdr > ref/mesh_cyl.hdr

mesh_cyl.hdr:	mesh.oct
	rad -v   cyl    mesh.rif

### End    mesh-cyl   tests

### Reference and test for    mesh view  cyl2 ###

test-mesh-cyl2:	ref/mesh_cyl2.hdr    mesh_cyl2.hdr
	$(RDU_PFILT)    mesh_cyl2.hdr | $(IMG_CMP) ref/mesh_cyl2.hdr -

ref/mesh_cyl2.hdr:
	make    mesh_cyl2.hdr
	$(RDU_PFILT)    mesh_cyl2.hdr > ref/mesh_cyl2.hdr

mesh_cyl2.hdr:	mesh.oct
	rad -v  cyl2    mesh.rif

### End    mesh-cyl2  tests

### Octree for  mirror ###

mirror.oct:	mirror.rif
	rad -v 0  mirror.rif

### Reference and test for  mirror view  fish ###

test-mirror-fish:	ref/mirror_fish.hdr  mirror_fish.hdr
	$(RDU_PFILT)  mirror_fish.hdr | $(IMG_CMP) ref/mirror_fish.hdr -

ref/mirror_fish.hdr:
	make  mirror_fish.hdr
	$(RDU_PFILT)  mirror_fish.hdr > ref/mirror_fish.hdr

mirror_fish.hdr:	mirror.oct
	rad -v  fish  mirror.rif

### End  mirror-fish  tests

### Octree for    mist ###

mist.oct:	mist.rif
	rad -v 0    mist.rif

### Reference and test for    mist view   def ###

test-mist-def:	ref/mist_def.hdr    mist_def.hdr
	$(RDU_PFILT)    mist_def.hdr | $(IMG_CMP) ref/mist_def.hdr -

ref/mist_def.hdr:
	make    mist_def.hdr
	$(RDU_PFILT)    mist_def.hdr > ref/mist_def.hdr

mist_def.hdr:	mist.oct
	rad -v   def    mist.rif

### End    mist-def   tests

### Octree for  prism1 ###

prism1.oct:	prism1.rif
	rad -v 0  prism1.rif

### Reference and test for  prism1 view  fish ###

test-prism1-fish:	ref/prism1_fish.hdr  prism1_fish.hdr
	$(RDU_PFILT)  prism1_fish.hdr | $(IMG_CMP) ref/prism1_fish.hdr -

ref/prism1_fish.hdr:
	make  prism1_fish.hdr
	$(RDU_PFILT)  prism1_fish.hdr > ref/prism1_fish.hdr

prism1_fish.hdr:	prism1.oct
	rad -v  fish  prism1.rif

### End  prism1-fish  tests

### Octree for  prism2 ###

prism2.oct:	prism2.rif
	rad -v 0  prism2.rif

### Reference and test for  prism2 view  fish ###

test-prism2-fish:	ref/prism2_fish.hdr  prism2_fish.hdr
	$(RDU_PFILT)  prism2_fish.hdr | $(IMG_CMP) ref/prism2_fish.hdr -

ref/prism2_fish.hdr:
	make  prism2_fish.hdr
	$(RDU_PFILT)  prism2_fish.hdr > ref/prism2_fish.hdr

prism2_fish.hdr:	prism2.oct
	rad -v  fish  prism2.rif

### End  prism2-fish  tests

### Octree for   tfunc ###

tfunc.oct:	tfunc.rif
	rad -v 0   tfunc.rif

### Reference and test for   tfunc view   def ###

test-tfunc-def:	ref/tfunc_def.hdr   tfunc_def.hdr
	$(RDU_PFILT)   tfunc_def.hdr | $(IMG_CMP) ref/tfunc_def.hdr -

ref/tfunc_def.hdr:
	make   tfunc_def.hdr
	$(RDU_PFILT)   tfunc_def.hdr > ref/tfunc_def.hdr

tfunc_def.hdr:	tfunc.oct
	rad -v   def   tfunc.rif

### End   tfunc-def   tests

### Reference and test for   tfunc view  fish ###

test-tfunc-fish:	ref/tfunc_fish.hdr   tfunc_fish.hdr
	$(RDU_PFILT)   tfunc_fish.hdr | $(IMG_CMP) ref/tfunc_fish.hdr -

ref/tfunc_fish.hdr:
	make   tfunc_fish.hdr
	$(RDU_PFILT)   tfunc_fish.hdr > ref/tfunc_fish.hdr

tfunc_fish.hdr:	tfunc.oct
	rad -v  fish   tfunc.rif

### End   tfunc-fish  tests

### Octree for   trans ###

trans.oct:	trans.rif
	rad -v 0   trans.rif

### Reference and test for   trans view   def ###

test-trans-def:	ref/trans_def.hdr   trans_def.hdr
	$(RDU_PFILT)   trans_def.hdr | $(IMG_CMP) ref/trans_def.hdr -

ref/trans_def.hdr:
	make   trans_def.hdr
	$(RDU_PFILT)   trans_def.hdr > ref/trans_def.hdr

trans_def.hdr:	trans.oct
	rad -v   def   trans.rif

### End   trans-def   tests

### Reference and test for   trans view  fish ###

test-trans-fish:	ref/trans_fish.hdr   trans_fish.hdr
	$(RDU_PFILT)   trans_fish.hdr | $(IMG_CMP) ref/trans_fish.hdr -

ref/trans_fish.hdr:
	make   trans_fish.hdr
	$(RDU_PFILT)   trans_fish.hdr > ref/trans_fish.hdr

trans_fish.hdr:	trans.oct
	rad -v  fish   trans.rif

### End   trans-fish  tests

### Octree for  trans2 ###

trans2.oct:	trans2.rif
	rad -v 0  trans2.rif

### Reference and test for  trans2 view   def ###

test-trans2-def:	ref/trans2_def.hdr  trans2_def.hdr
	$(RDU_PFILT)  trans2_def.hdr | $(IMG_CMP) ref/trans2_def.hdr -

ref/trans2_def.hdr:
	make  trans2_def.hdr
	$(RDU_PFILT)  trans2_def.hdr > ref/trans2_def.hdr

trans2_def.hdr:	trans2.oct
	rad -v   def  trans2.rif

### End  trans2-def   tests

### Reference and test for  trans2 view  fish ###

test-trans2-fish:	ref/trans2_fish.hdr  trans2_fish.hdr
	$(RDU_PFILT)  trans2_fish.hdr | $(IMG_CMP) ref/trans2_fish.hdr -

ref/trans2_fish.hdr:
	make  trans2_fish.hdr
	$(RDU_PFILT)  trans2_fish.hdr > ref/trans2_fish.hdr

trans2_fish.hdr:	trans2.oct
	rad -v  fish  trans2.rif

### End  trans2-fish  tests

### Reference and test for  trans2 view   win ###

test-trans2-win:	ref/trans2_win.hdr  trans2_win.hdr
	$(RDU_PFILT)  trans2_win.hdr | $(IMG_CMP) ref/trans2_win.hdr -

ref/trans2_win.hdr:
	make  trans2_win.hdr
	$(RDU_PFILT)  trans2_win.hdr > ref/trans2_win.hdr

trans2_win.hdr:	trans2.oct
	rad -v   win  trans2.rif

### End  trans2-win   tests
