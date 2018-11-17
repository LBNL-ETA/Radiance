# RCSid $Id: Makefile,v 1.1 2018/11/17 22:09:12 greg Exp $
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
	rm -f *.oct *.amb *.dat *_*.hdr

### Reference and test for dielectric view def ###

test-dielectric-def:	ref/dielectric_def.hdr dielectric_def.hdr
	$(RDU_PFILT) dielectric_def.hdr | $(IMG_CMP) ref/dielectric_def.hdr -

ref/dielectric_def.hdr:
	make dielectric_def.hdr
	$(RDU_PFILT) dielectric_def.hdr > ref/dielectric_def.hdr

dielectric_def.hdr:	dielectric.rif
	rad -v def dielectric.rif

### End dielectric-def tests

### Reference and test for dielectric view fish ###

test-dielectric-fish:	ref/dielectric_fish.hdr dielectric_fish.hdr
	$(RDU_PFILT) dielectric_fish.hdr | $(IMG_CMP) ref/dielectric_fish.hdr -

ref/dielectric_fish.hdr:
	make dielectric_fish.hdr
	$(RDU_PFILT) dielectric_fish.hdr > ref/dielectric_fish.hdr

dielectric_fish.hdr:	dielectric.rif
	rad -v fish dielectric.rif

### End dielectric-fish tests

### Reference and test for glass view def ###

test-glass-def:	ref/glass_def.hdr glass_def.hdr
	$(RDU_PFILT) glass_def.hdr | $(IMG_CMP) ref/glass_def.hdr -

ref/glass_def.hdr:
	make -s glass_def.hdr
	$(RDU_PFILT) glass_def.hdr > ref/glass_def.hdr

glass_def.hdr:	glass.rif
	rad -v def glass.rif

### End glass-def tests

### Reference and test for glass view fish ###

test-glass-fish:	ref/glass_fish.hdr glass_fish.hdr
	$(RDU_PFILT) glass_fish.hdr | $(IMG_CMP) ref/glass_fish.hdr -

ref/glass_fish.hdr:
	make -s glass_fish.hdr
	$(RDU_PFILT) glass_fish.hdr > ref/glass_fish.hdr

glass_fish.hdr:	glass.rif
	rad -v fish glass.rif

### End glass-fish tests

### Reference and test for glass view up ###

test-glass-up:	ref/glass_up.hdr glass_up.hdr
	$(RDU_PFILT) glass_up.hdr | $(IMG_CMP) ref/glass_up.hdr -

ref/glass_up.hdr:
	make -s glass_up.hdr
	$(RDU_PFILT) glass_up.hdr > ref/glass_up.hdr

glass_up.hdr:	glass.rif
	rad -v up glass.rif

### End glass-up tests

### Reference and test for inst view def ###

test-inst-def:	ref/inst_def.hdr inst_def.hdr
	$(RDU_PFILT) inst_def.hdr | $(IMG_CMP) ref/inst_def.hdr -

ref/inst_def.hdr:
	make -s inst_def.hdr
	$(RDU_PFILT) inst_def.hdr > ref/inst_def.hdr

inst_def.hdr:	inst.rif
	rad -v def inst.rif

### End inst-def tests

### Reference and test for inst view fish ###

test-inst-fish:	ref/inst_fish.hdr inst_fish.hdr
	$(RDU_PFILT) inst_fish.hdr | $(IMG_CMP) ref/inst_fish.hdr -

ref/inst_fish.hdr:
	make -s inst_fish.hdr
	$(RDU_PFILT) inst_fish.hdr > ref/inst_fish.hdr

inst_fish.hdr:	inst.rif
	rad -v fish inst.rif

### End inst-fish tests

### Reference and test for mesh view def ###

test-mesh-def:	ref/mesh_def.hdr mesh_def.hdr
	$(RDU_PFILT) mesh_def.hdr | $(IMG_CMP) ref/mesh_def.hdr -

ref/mesh_def.hdr:
	make -s mesh_def.hdr
	$(RDU_PFILT) mesh_def.hdr > ref/mesh_def.hdr

mesh_def.hdr:	mesh.rif
	rad -v def mesh.rif

### End mesh-def tests

### Reference and test for mesh view cyl ###

test-mesh-cyl:	ref/mesh_cyl.hdr mesh_cyl.hdr
	$(RDU_PFILT) mesh_cyl.hdr | $(IMG_CMP) ref/mesh_cyl.hdr -

ref/mesh_cyl.hdr:
	make -s mesh_cyl.hdr
	$(RDU_PFILT) mesh_cyl.hdr > ref/mesh_cyl.hdr

mesh_cyl.hdr:	mesh.rif
	rad -v cyl mesh.rif

### End mesh-cyl tests

### Reference and test for mesh view cyl2 ###

test-mesh-cyl2:	ref/mesh_cyl2.hdr mesh_cyl2.hdr
	$(RDU_PFILT) mesh_cyl2.hdr | $(IMG_CMP) ref/mesh_cyl2.hdr -

ref/mesh_cyl2.hdr:
	make -s mesh_cyl2.hdr
	$(RDU_PFILT) mesh_cyl2.hdr > ref/mesh_cyl2.hdr

mesh_cyl2.hdr:	mesh.rif
	rad -v cyl2 mesh.rif

### End mesh-cyl2 tests

### Reference and test for mirror view fish ###

test-mirror-fish:	ref/mirror_fish.hdr mirror_fish.hdr
	$(RDU_PFILT) mirror_fish.hdr | $(IMG_CMP) ref/mirror_fish.hdr -

ref/mirror_fish.hdr:
	make -s mirror_fish.hdr
	$(RDU_PFILT) mirror_fish.hdr > ref/mirror_fish.hdr

mirror_fish.hdr:	mirror.rif
	rad -v fish mirror.rif

### End mirror-fish tests

### Reference and test for mist view def ###

test-mist-def:	ref/mist_def.hdr mist_def.hdr
	$(RDU_PFILT) mist_def.hdr | $(IMG_CMP) ref/mist_def.hdr -

ref/mist_def.hdr:
	make -s mist_def.hdr
	$(RDU_PFILT) mist_def.hdr > ref/mist_def.hdr

mist_def.hdr:	mist.rif
	rad -v def mist.rif

### End mist-def tests

### Reference and test for prism1 view fish ###

test-prism1-fish:	ref/prism1_fish.hdr prism1_fish.hdr
	$(RDU_PFILT) prism1_fish.hdr | $(IMG_CMP) ref/prism1_fish.hdr -

ref/prism1_fish.hdr:
	make -s prism1_fish.hdr
	$(RDU_PFILT) prism1_fish.hdr > ref/prism1_fish.hdr

prism1_fish.hdr:	prism1.rif
	rad -v fish prism1.rif

### End prism1-fish tests

### Reference and test for prism2 view fish ###

test-prism2-fish:	ref/prism2_fish.hdr prism2_fish.hdr
	$(RDU_PFILT) prism2_fish.hdr | $(IMG_CMP) ref/prism2_fish.hdr -

ref/prism2_fish.hdr:
	make -s prism2_fish.hdr
	$(RDU_PFILT) prism2_fish.hdr > ref/prism2_fish.hdr

prism2_fish.hdr:	prism2.rif
	rad -v fish prism2.rif

### End prism2-fish tests

### Reference and test for tfunc view def ###

test-tfunc-def:	ref/tfunc_def.hdr tfunc_def.hdr
	$(RDU_PFILT) tfunc_def.hdr | $(IMG_CMP) ref/tfunc_def.hdr -

ref/tfunc_def.hdr:
	make -s tfunc_def.hdr
	$(RDU_PFILT) tfunc_def.hdr > ref/tfunc_def.hdr

tfunc_def.hdr:	tfunc.rif
	rad -v def tfunc.rif

### End tfunc-def tests

### Reference and test for tfunc view fish ###

test-tfunc-fish:	ref/tfunc_fish.hdr tfunc_fish.hdr
	$(RDU_PFILT) tfunc_fish.hdr | $(IMG_CMP) ref/tfunc_fish.hdr -

ref/tfunc_fish.hdr:
	make -s tfunc_fish.hdr
	$(RDU_PFILT) tfunc_fish.hdr > ref/tfunc_fish.hdr

tfunc_fish.hdr:	tfunc.rif
	rad -v fish tfunc.rif

### End tfunc-fish tests

### Reference and test for trans view def ###

test-trans-def:	ref/trans_def.hdr trans_def.hdr
	$(RDU_PFILT) trans_def.hdr | $(IMG_CMP) ref/trans_def.hdr -

ref/trans_def.hdr:
	make -s trans_def.hdr
	$(RDU_PFILT) trans_def.hdr > ref/trans_def.hdr

trans_def.hdr:	trans.rif
	rad -v def trans.rif

### End trans-def tests

### Reference and test for trans view fish ###

test-trans-fish:	ref/trans_fish.hdr trans_fish.hdr
	$(RDU_PFILT) trans_fish.hdr | $(IMG_CMP) ref/trans_fish.hdr -

ref/trans_fish.hdr:
	make -s trans_fish.hdr
	$(RDU_PFILT) trans_fish.hdr > ref/trans_fish.hdr

trans_fish.hdr:	trans.rif
	rad -v fish trans.rif

### End trans-fish tests

### Reference and test for trans2 view def ###

test-trans2-def:	ref/trans2_def.hdr trans2_def.hdr
	$(RDU_PFILT) trans2_def.hdr | $(IMG_CMP) ref/trans2_def.hdr -

ref/trans2_def.hdr:
	make -s trans2_def.hdr
	$(RDU_PFILT) trans2_def.hdr > ref/trans2_def.hdr

trans2_def.hdr:	trans2.rif
	rad -v def trans2.rif

### End trans2-def tests

### Reference and test for trans2 view fish ###

test-trans2-fish:	ref/trans2_fish.hdr trans2_fish.hdr
	$(RDU_PFILT) trans2_fish.hdr | $(IMG_CMP) ref/trans2_fish.hdr -

ref/trans2_fish.hdr:
	make -s trans2_fish.hdr
	$(RDU_PFILT) trans2_fish.hdr > ref/trans2_fish.hdr

trans2_fish.hdr:	trans2.rif
	rad -v fish trans2.rif

### End trans2-fish tests

### Reference and test for trans2 view win ###

test-trans2-win:	ref/trans2_win.hdr trans2_win.hdr
	$(RDU_PFILT) trans2_win.hdr | $(IMG_CMP) ref/trans2_win.hdr -

ref/trans2_win.hdr:
	make -s trans2_win.hdr
	$(RDU_PFILT) trans2_win.hdr > ref/trans2_win.hdr

trans2_win.hdr:	trans2.rif
	rad -v win trans2.rif

### End trans2-win tests
