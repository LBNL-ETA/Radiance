#!/bin/csh -f
# RCSid: $Id: objpict.csh,v 2.3 2004/01/29 22:45:30 greg Exp $
#
# Make a nice multi-view picture of an object
# Command line arguments contain materials and object files
#
set objdir=/usr/local/lib/ray/lib
set tmpdir=/usr/tmp
set xres=250
set yres=250
set rpict="rpict -av .2 .2 .2 -x $xres -y $yres"
set inprad=$tmpdir/op$$.rad
set octree=$tmpdir/op$$.oct
set pict1=$tmpdir/op$$a.pic
set pict2=$tmpdir/op$$b.pic
set pict3=$tmpdir/op$$c.pic
set pict4=$tmpdir/op$$d.pic
set tmpfiles="$inprad $octree $pict1 $pict2 $pict3 $pict4"
onintr quit
if ( $#argv ) then
	cat $* > $inprad
else
	cat > $inprad
endif
set dims=`getbbox -h $inprad`
set siz=`rcalc -n -e 'max(a,b):if(a-b,a,b);$1='"max($dims[2]-$dims[1],max($dims[4]-$dims[3],$dims[6]-$dims[5]))"`
set vw1="-vtl -vp 2 .5 .5 -vd -1 0 0 -vh 1 -vv 1"
set vw2="-vtl -vp .5 2 .5 -vd 0 -1 0 -vh 1 -vv 1"
set vw3="-vtl -vp .5 .5 2 -vd 0 0 -1 -vu -1 0 0 -vh 1 -vv 1"
set vw4="-vp 3 3 3 -vd -1 -1 -1 -vh 20 -vv 20"

xform -t `ev "-($dims[1]+$dims[2])/2" "-($dims[3]+$dims[4])/2" "-($dims[5]+$dims[6])/2"` \
		-s `ev 1/$siz` -t .5 .5 .5 $inprad \
	| oconv $objdir/testroom - > $octree
$rpict $vw1 $octree > $pict1
$rpict $vw2 $octree > $pict2
$rpict $vw3 $octree > $pict3
$rpict $vw4 $octree > $pict4
pcompos $pict3 0 $yres $pict4 $xres $yres $pict1 0 0 $pict2 $xres 0

quit:
rm -f $tmpfiles
exit 0
