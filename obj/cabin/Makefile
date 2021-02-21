#
# Makefile for Cabin Scene
#

VIEW= -vf vf/plan

SCENE= summercabin

AMB= -av .1 .1 .1

view:	$(SCENE).oct
	rvu $(VIEW) $(AMB) $(SCENE).oct

sholo:	living.oct bed.oct bath.oct summer.hdk
	rholo -n 1 -o ogl summer.hdk

nholo:	living.oct bed.oct bath.oct night.hdk
	rholo -n 1 -o ogl night.hdk

clean:
	set nonomatch; rm -f *.hdk *.oct *.amb

summer.hdk:	summer.hif summercabin.oct
	rholo -f summer.hdk summer.hif

night.hdk:	night.hif nightcabin.oct
	rholo -f night.hdk night.hif

summercabin.oct:	cabin.oct summerday landscape lights.off \
daywindows pattmats
	oconv -f -r 8192 -i cabin.oct summerday landscape lights.off \
daywindows > summercabin.oct

wintercabin.oct:	cabin.oct winterday landscape daywindows \
lights.off pattmats
	oconv -f -r 8192 -i cabin.oct winterday landscape lights.off \
daywindows > wintercabin.oct

living.oct:	testmats cabin lights.off lvfurn.rad
	oconv -f testmats cabin lights.off lvfurn.rad > living.oct

bed.oct:	testmats cabin lights.off lvfurn.rad
	oconv -f testmats cabin lights.off bdfurn.rad > bed.oct

bath.oct:	testmats cabin bathroom lights.off
	oconv -f testmats cabin lights.off bathroom > bath.oct

nightcabin.oct:	cabin.oct lights pattmats
	oconv -f -r 8192 -i cabin.oct lights > nightcabin.oct

cabin.oct:	cabin bathroom furniture winpanes.rad mirrors.rad
	oconv -b -100 -100 -100 225 -r 8192 \
pattmats cabin bathroom furniture winpanes.rad mirrors.rad > cabin.oct

cabin.oct:	chair.norm coatrack.norm door.norm dresser.norm sofa.norm \
winframe.norm

living.oct:	chair.norm coatrack.norm door.norm sofa.norm winframe.norm

bed.oct:	chair.norm door.norm dresser.norm winframe.norm

nightcabin.oct:	flood.rad

summercabin.oct wintercabin.oct:	flood.shroud
