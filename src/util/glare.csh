#!/bin/csh -f
# RCSid: $Id: glare.csh,v 2.4 2003/02/22 02:07:30 greg Exp $
#
# Interactive program for calculating glare values
#
set fgargs=(-v)
set nofile="none"
set octree=$nofile
set picture=$nofile
set glarefile=$nofile
if ($#argv >= 1) then
	set glarefile=$argv[1]
	if ($#argv >= 2) then
		set picture=$argv[2]
		if ($#argv >= 3) then
			set octree=$argv[3]
		endif
	endif
endif
set rtargs=
set view=
set threshold=300
set maxangle=0
set stepangle=10

set gndxfile="glr$$.ndx"
set plotfile="glr$$.plt"
set tempfiles=($gndxfile $plotfile)

alias readvar 'echo -n Enter \!:1 "[$\!:1]: ";set ans="$<";if("$ans" != "")set \!:1="$ans"'

onintr quit

cat <<_EOF_
This script provides a convenient interface to the glare calculation
and analysis programs.  It works on a single output file from findglare,
which you must either create here or have created already before.
_EOF_
while ( $glarefile == $nofile )
	echo ""
	echo "Please specify a glare file, or use ^C to exit."
	readvar glarefile
end
if ( ! -f $glarefile ) then
	echo ""
	echo -n "The glare file '$glarefile' doesn't exist -- shall we create it? "
	if ("$<" =~ [nN]*) exit 0
	cat <<_EOF_

A glare file can be created from a wide-angle Radiance picture (preferably
a fisheye view), an octree, or both.  Computing glare source locations
from a finished image is faster, but using an octree is more reliable
since all of the scene information is there.  Whenever it is possible,
you should use both a picture and an octree to compute the glare file.

You must also have a viewpoint and direction for the glare calculation.
If the viewpoint is not the same as the picture, an error will result.
If no view is specified, the view parameters are taken from the picture,
so the view specification is optional unless you are starting from an
octree.

_EOF_
	if ( $picture == $nofile ) then
		readvar picture
	endif
	if ( $octree == $nofile ) then
		readvar octree
	endif
	readvar view
	if ( $picture == $nofile && $octree == $nofile ) then
		echo "You must specify a picture or an octree"
		exit 1
	endif
	if ( $picture == $nofile && "$view" == "" ) then
		echo "You must give a view if there is no picture"
		exit 1
	endif
	if ( $picture != $nofile ) then
		if ( ! -r $picture ) then
			echo "Cannot read $picture"
			exit 1
		endif
		set fgargs=($fgargs -p $picture)
		if ( "$view" != "" ) then
			set view=(-vf $picture $view)
		endif
	endif
	set fgargs=($fgargs $view)
	if ( $octree != $nofile ) then
		if ( ! -r $octree ) then
			echo "Cannot read $octree"
			exit 1
		endif
		cat <<_EOF_

With an octree, you should give the same options for -av, -ab and
so forth as are used to render the scene.  Please enter them below.

_EOF_
		if ($picture != $nofile) then
			echo "These are the parameters from the picture:"
			getinfo $picture
		endif
		readvar rtargs
		set fgargs=($fgargs $rtargs $octree)
	endif
	echo ""
	echo -n "Do you expect glare sources in the scene to be small? "
	if ("$<" =~ [yY]) then
		echo -n "Very small? "
		if ("$<" =~ [yY]) then
			set fgargs=(-r 400 $fgargs)
		else
			set fgargs=(-r 250 $fgargs)
		endif
		echo -n "Are these sources relatively bright? "
		if ("$<" =~ [nN]*) then
			set fgargs=(-c $fgargs)
		endif
	endif
	cat <<_EOF_

Glare sources are determined by a threshold level.  Any part of the
view that is above this threshold is considered to be part of a glare
source.  If you do not choose a threshold value, it will be set auto-
matically to 7 times the average field luminance (in candelas/sq.meter).

_EOF_
	echo -n "Do you wish to set the glare threshold manually? "
	if ("$<" =~ [yY]) then
		readvar threshold
		set fgargs=(-t $threshold $fgargs)
	endif
	cat <<_EOF_

It is often desirable to know how glare changes as a function of
viewing direction.  If you like, you can compute glare for up
to 180 degrees to the left and right of the view center.  Enter
the maximum angle that you would like to compute below.

_EOF_
	readvar maxangle
	if ( $maxangle >= $stepangle ) then
		set fgargs=(-ga $stepangle-$maxangle\:$stepangle $fgargs)
	endif
	echo ""
	echo "Starting calculation..."
	echo findglare $fgargs
	findglare $fgargs > $glarefile
	if ($status) then
		echo "Uh-oh.  Something went wrong with findglare\!"
		rm $glarefile
		exit 1
	endif
endif
if ($?DISPLAY && $picture != $nofile) then
	echo ""
	echo "Displaying glare sources in '$picture'..."
	xglaresrc $picture $glarefile
endif
set ndxnam=("Guth Visual Comfort Probability" "Guth Disability Glare Ratio" "CIE Glare Index" "BRS Glare Index" "Unified Glare Rating" "Daylight Glare Index")
set ndxcom=("glarendx -t guth_vcp" "glarendx -t guth_dgr" "glarendx -t cie_cgi" "glarendx -t brs_gi" "glarendx -t ugr" "glarendx -t dgi")
set bcdlvl=(50 124 18.5 18.5 18.5 21.7)
while ( 1 )
	echo ""
	echo "	0.  Quit"
	set i=1
	while ($i <= $#ndxnam)
		echo "	$i.  $ndxnam[$i]"
		@ i++
	end
	echo ""
	set choice=0
	readvar choice
	if ($choice < 1 || $choice > $#ndxcom) goto quit
	$ndxcom[$choice] $glarefile > $gndxfile
	echo ""
	echo $ndxnam[$choice]
	echo "Angle		Value"
	getinfo - < $gndxfile
	echo ""
	set save_file=$nofile
	readvar save_file
	if ($save_file != $nofile) then
		cp -i $gndxfile $save_file
	endif
	if ( `getinfo - < $gndxfile | wc -l` > 4 ) then
		echo ""
		echo -n "Do you want to plot this data? "
		if ( "$<" =~ [yY]* ) then
			set subtitle="$view"
			if ($picture != $nofile) then
				set subtitle="$picture $subtitle"
			else if ($octree != $nofile) then
				set subtitle="$subtitle $octree"
			endif
			if ( "$subtitle" == "" ) then
				set subtitle="`getinfo < $glarefile | grep findglare`"
			endif
			cat <<_EOF_ > $plotfile
include=line.plt
include=polar.plt
title="$ndxnam[$choice]"
subtitle="$subtitle"
Bdata=
	0	$bcdlvl[$choice]
	360	$bcdlvl[$choice]
;
Blabel="BCD level"
Adata=
_EOF_
			getinfo - < $gndxfile >> $plotfile
			igraph $plotfile
		endif
	endif
end
quit:
rm -f $tempfiles
exit 0
