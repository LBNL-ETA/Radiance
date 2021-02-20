#
# Makefile for our little office scene
#

DEV=x11
VW=outside
VIEW=-vf $(VW).vp
SCENE=plain
AMB=-av 4 4 4
DIR=-dr 2 -dp 64 -dj .5
RPOPTS=-t 1800 -ps 1

view:	$(SCENE).oct
	rvu $(VIEW) -o $(DEV) $(AMB) $(DIR) $(SCENE).oct

picture:	$(SCENE)_$(VW).hdr

$(SCENE)_$(VW).hdr:	$(SCENE).oct
	rpict $(VIEW) $(AMB) $(DIR) $(RPOPTS) \
$(SCENE).oct > $(SCENE)_$(VW).unf
	mv $(SCENE)_$(VW).unf $(SCENE)_$(VW).hdr

recover:	$(SCENE)_$(VW).unf
	mv $(SCENE)_$(VW).unf $(SCENE)_$(VW).rvr
	rpict $(AMB) $(DIR) $(RPOPTS) -r $(SCENE)_$(VW).rvr \
$(SCENE).oct > $(SCENE)_$(VW).unf
	mv $(SCENE)_$(VW).unf $(SCENE)_$(VW).hdr

plainfurn.oct:	plain.oct furniture.rad desk.rad draft_ang.rad draft_flat.rad \
orange_chair.rad stool.rad wastebasket.rad
	oconv -f -i plain.oct furniture.rad > plainfurn.oct

plain.oct:	plainoff.rad environ.rad
	oconv plainoff.rad environ.rad > plain.oct

illumfurn.oct:	illum.oct furniture.rad desk.rad draft_ang.rad draft_flat.rad \
orange_chair.rad stool.rad wastebasket.rad
	oconv -f -i illum.oct furniture.rad > illumfurn.oct

illum.oct:	illumoff.rad environ.rad
	oconv illumoff.rad environ.rad > illum.oct

illumoff.rad:	plainoff.rad plain.oct
	mkillum -dr 2 -ab 1 -ad 128 -as 64 -av .1 .1 .1 plain.oct \
< plainoff.rad > illumoff.rad
