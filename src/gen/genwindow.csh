#!/bin/csh -f
# RCSid: $Id: genwindow.csh,v 2.2 2003/02/22 02:07:23 greg Exp $
#
#  Generate source description for window with venetian blinds.
#
#	1/11/88		Greg Ward
#
#  Usage: genwindow shellfile
#
#  Takes as input the following shell variables from file $1:
#
#	worient		Window orientation (degrees from south)
#	wwidth		Window width
#	wheight		Window height
#	wtrans		Window transmission
#	bdepth		Blind depth
#	bspac		Blind vertical spacing
#	gap		Gap between blinds and window
#	brcurv		Blind curvature radius (negative upward, zero none)
#	bangle		Blind inclination (degrees from horizontal, up and out)
#	material	Blind material type (metal or plastic)
#	ucolor		Blind up side color ("red green blue")
#	uspec		Blind up side specularity (0 to 1)
#	urough		Blind up side roughness (0 to .2)
#	dcolor		Blind down side color (opt)
#	dspec		Blind down side specularity (opt)
#	drough		Blind down side roughness (opt)
#	lat		Site latitude (degrees)
#	long		Site longitude (degrees)
#	mer		Site standard meridian (degrees)
#	hour		Hour (standard 24hr time)
#	day		Day (of month)
#	month		Month (of year)
#	grefl		Ground plane reflectance
#	sky		Sky conditions (sunny, clear, cloudy)
#	nsources	Number of sources to divide window
#
#  Creates the following output files:
#
#	stdout		Window description for Radiance
#	$1.d$$		Window output distribution
#

if ( $#argv != 1 ) then
	echo Usage: $0 input
	exit 1
endif

#	Set input and output files
set input = $1
set distrib = $1.d$$
set tmpdir = /tmp
set descrip = $tmpdir/gw$$.des
set distoct = $tmpdir/gw$$.oct
set skydesc = $tmpdir/gw$$.sky
set remove = ( $distoct $skydesc $descrip )
set removerr = ( $remove $distrib )
onintr error

#	Set default values
set worient = 0.
set wtrans = ( .96 .96 .96 )
set bdepth = 0.025
set bspac = 0.0155
set gap = 0.
set brcurv = 0.
set bangle = 0.
set material = plastic
set ucolor = ( 0.5 0.5 0.5 ) 
set uspec = 0.
set urough = 0.05
set lat = 37.8
set long = 122
set mer = 120
set hour = 12
set day = 21
set month = 10
set grefl = .2
set sky = sunny
set nsources = 6

#	Get input
source $input
if ( $status ) goto error

#	Create window
cat > $descrip <<_EOF_

void glass clear_glass
0
0
3 $wtrans

clear_glass polygon window
0
0
12
	0	0	0
	0	0	$wheight
	$wwidth	0	$wheight
	$wwidth	0	0

void $material blind_upmat
0
0
5 $ucolor $uspec $urough
_EOF_
if ( $status ) goto error

#	Blinds
genblinds blind_upmat blind_up $bdepth $wwidth $wheight \
		`ev "ceil($wheight/$bspac)"` $bangle -r $brcurv \
	| xform -t $gap -$wwidth 1e-4 -rz 90 >> $descrip
if ( $status ) goto error
if ( $?dcolor ) then
	cat >> $descrip <<_EOF_

	void $material blind_dnmat
	0
	0
	5 $dcolor $dspec $drough
_EOF_
	genblinds blind_dnmat blind_down $bdepth $wwidth $wheight \
			`ev "ceil($wheight/$bspac)"` $bangle -r $brcurv \
		| xform -t $gap -$wwidth 0 -rz 90 >> $descrip
	if ( $status ) goto error
endif

#	Make sky
switch ($sky)
case sun*:
	set skysw = +s
	breaksw
case clear:
	set skysw = -s
	breaksw
case cloud*:
	set skysw = -c
	breaksw
endsw
dosky:
gensky $month $day $hour $skysw -a $lat -o $long -m $mer -g $grefl \
	| xform -rz `ev "-($worient)"` > $skydesc
if ( $skysw == +s ) then
	if ( `sed -n 13p $skydesc | rcalc -e '$1=if($3,-1,1)'` < 0 ) then
		set skysw = -s
		goto dosky
	endif
endif
cat >> $skydesc <<_EOF_

skyfunc glow skyglow
0
0
4 1 1 1 0

skyglow source sky
0
0
4 0 -1 0 180
_EOF_

#	Make distribution
oconv $skydesc $descrip > $distoct
if ( $status ) goto error
echo 2 5 85 9 0 340 18 > $distrib
makedist -h -d -x1 0 1 0 -x2 0 0 1 -x3 1 0 0 -alpha 5-85:10 -beta 0-340:20 \
		-tw $bspac -th $bspac -td `ev "2*($bdepth+$gap)"` \
		-tc `ev "$wwidth/2" "$bdepth+$gap" "$wheight/2"` \
		-ab 2 -aa .1 -ad 64 -as 64 \
		-x 16 -y 16 $distoct >> $distrib
if ( $status ) goto error

set wsgrid = (`ev "floor(sqrt($nsources*$wwidth/$wheight)+.5)" "floor(sqrt($nsources*$wheight/$wwidth)+.5)"`)
@ nsources = $wsgrid[1] * $wsgrid[2]
#	Print header
cat <<_EOF_
#
#  Window with venetian blinds
#  Created from $input `date`
#
#	Window orientation (degrees from south):	$worient
#	Window height:					$wheight
#	Window width:					$wwidth
#	Window transmission:				$wtrans
#	Blind depth:					$bdepth
#	Blind spacing:					$bspac
#	Gap to window:					$gap
#	Blind curvature radius:				$brcurv
#	Blind inclination (degrees altitude):		$bangle
#	Blind material:					$material
#	Up side color:					$ucolor
#	Up side specularity:				$uspec
#	Up side roughness:				$urough
_EOF_
if ( $?dcolor ) then
	cat <<_EOF_
#	Down side color:				$dcolor
#	Down side specularity:				$dspec
#	Down side roughness:				$drough
_EOF_
else
	echo \#\	Down side same as up
endif
cat <<_EOF_
#	Latitude (degrees):				$lat
#	Longitude (degrees):				$long
#	Standard Meridian (degrees):			$mer
#	Month Day Hour:					$month $day $hour
#	Ground plane reflectance:			$grefl
#	Sky condition:					$sky
#	Number of window sources:			$nsources
#
_EOF_

#	Send sources
xform -e -rz $worient <<_EOF_

void brightdata wdistrib
10 noop $distrib source.cal src_theta src_phi -rx 90 -ry -90 -mx
0
0

wdistrib illum willum
0
0
3 1 1 1

!gensurf willum wsource "$wwidth*t" $bdepth "$wheight*s" $wsgrid[2] $wsgrid[1]
_EOF_
if ( $status ) goto error

#	Send window
xform -rz $worient $descrip
if ( $status ) goto error

#	All done, print and exit
rm -f $remove
exit 0

#	Error exit
error:
	rm -f $removerr
	exit 1
