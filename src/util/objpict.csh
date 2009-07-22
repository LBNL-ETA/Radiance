#!/bin/csh -f
# RCSid: $Id: objpict.csh,v 2.6 2009/07/22 17:35:06 greg Exp $
#
# Make a nice multi-view picture of an object
# Command line arguments contain materials and object files
#
set tmpdir=/tmp/objv$$
set xres=250
set yres=250
set rpict="rpict -av .2 .2 .2 -x $xres -y $yres"
set inprad=$tmpdir/op.rad
set testroom=$tmpdir/testroom.rad
set octree=$tmpdir/op.oct
set pict1=$tmpdir/opa.hdr
set pict2=$tmpdir/opb.hdr
set pict3=$tmpdir/opc.hdr
set pict4=$tmpdir/opd.hdr
onintr quit
if ( $#argv ) then
	cat $* > $inprad
else
	cat > $inprad
endif
cat > $testroom << '_EOF_'
void plastic wall_mat 0 0 5 .681 .543 .686 0 .2
void light bright 0 0 3 3000 3000 3000
bright sphere lamp0 0 0 4 4 4 -4 .1
bright sphere lamp1 0 0 4 4 0 4 .1
bright sphere lamp2 0 0 4 0 4 4 .1
wall_mat polygon box.1540 0 0 12
                  5                 -5                 -5
                  5                 -5                  5
                 -5                 -5                  5
                 -5                 -5                 -5
wall_mat polygon box.4620 0 0 12
                 -5                 -5                  5
                 -5                  5                  5
                 -5                  5                 -5
                 -5                 -5                 -5
wall_mat polygon box.2310 0 0 12
                 -5                  5                 -5
                  5                  5                 -5
                  5                 -5                 -5
                 -5                 -5                 -5
wall_mat polygon box.3267 0 0 12
                  5                  5                 -5
                 -5                  5                 -5
                 -5                  5                  5
                  5                  5                  5
wall_mat polygon box.5137 0 0 12
                  5                 -5                  5
                  5                 -5                 -5
                  5                  5                 -5
                  5                  5                  5
wall_mat polygon box.6457 0 0 12
                 -5                  5                  5
                 -5                 -5                  5
                  5                 -5                  5
                  5                  5                  5
'_EOF_'
set dims=`getbbox -h $inprad`
set siz=`rcalc -n -e 'max(a,b):if(a-b,a,b);$1='"max($dims[2]-$dims[1],max($dims[4]-$dims[3],$dims[6]-$dims[5]))"`
set vw1="-vtl -vp 2 .5 .5 -vd -1 0 0 -vh 1 -vv 1"
set vw2="-vtl -vp .5 2 .5 -vd 0 -1 0 -vh 1 -vv 1"
set vw3="-vtl -vp .5 .5 2 -vd 0 0 -1 -vu -1 0 0 -vh 1 -vv 1"
set vw4="-vp 3 3 3 -vd -1 -1 -1 -vh 20 -vv 20"

xform -t `ev "-($dims[1]+$dims[2])/2" "-($dims[3]+$dims[4])/2" "-($dims[5]+$dims[6])/2"` \
		-s `ev 1/$siz` -t .5 .5 .5 $inprad \
	| oconv $testroom - > $octree
$rpict $vw1 $octree > $pict1
$rpict $vw2 $octree > $pict2
$rpict $vw3 $octree > $pict3
$rpict $vw4 $octree > $pict4
pcompos $pict3 0 $yres $pict4 $xres $yres $pict1 0 0 $pict2 $xres 0

quit:
rm -r $tmpdir
exit 0
