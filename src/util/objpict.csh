#!/bin/csh -f
# RCSid: $Id: objpict.csh,v 2.2 2003/02/22 02:07:30 greg Exp $
#
# Make a nice multi-view picture of an object
# Standard input is description made to fit in unit cube in positive quadrant
#
set objdir=/usr/local/lib/ray/lib
set tmpdir=/usr/tmp
set xres=250
set yres=250
set rpict="rpict -av .2 .2 .2 -x $xres -y $yres $*"
set octree=$tmpdir/op$$.oct
set pict1=$tmpdir/op$$a.pic
set pict2=$tmpdir/op$$b.pic
set pict3=$tmpdir/op$$c.pic
set pict4=$tmpdir/op$$d.pic
set tmpfiles="$octree $pict1 $pict2 $pict3 $pict4"
onintr quit
set vw1="-vtl -vp 2 .5 .5 -vd -1 0 0 -vh 1 -vv 1"
set vw2="-vtl -vp .5 2 .5 -vd 0 -1 0 -vh 1 -vv 1"
set vw3="-vtl -vp .5 .5 2 -vd 0 0 -1 -vu -1 0 0 -vh 1 -vv 1"
set vw4="-vp 3 3 3 -vd -1 -1 -1 -vh 20 -vv 20"

oconv $objdir/testroom - > $octree
$rpict $vw1 $octree > $pict1
$rpict $vw2 $octree > $pict2
$rpict $vw3 $octree > $pict3
$rpict $vw4 $octree > $pict4
pcompos $pict3 0 $yres $pict4 $xres $yres $pict1 0 0 $pict2 $xres 0

quit:
rm -f $tmpfiles
exit 0
