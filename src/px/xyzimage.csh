#!/bin/csh -f
# RCSid: $Id: xyzimage.csh,v 2.5 2008/08/25 04:50:32 greg Exp $
#
# Display one or more CIE XYZE pictures using ximage
#
set popt=""
if ( $?DISPLAY_PRIMARIES ) then
	set popt="-p $DISPLAY_PRIMARIES"
endif
set xiargs=""
set i=1
set firstarg=0
while ( $i <= $#argv && ! $firstarg )
	@ i1= $i + 1
	switch ( $argv[$i] )
	case -ge*:
	case -di*:
	case -g:
	case -c:
	case -e:
		if ( $i1 > $#argv ) goto notenough
		set xiargs=($xiargs $argv[$i] $argv[$i1])
		@ i1++
		breaksw
	case -d:
	case -b:
	case -m:
	case -f:
	case -s:
	case =*:
	case -o*:
		set xiargs=($xiargs $argv[$i])
		breaksw
	case -p:
		set popt="-p"
		while ( $i1 <= $i + 8 )
			if ( $i1 > $#argv ) goto notenough
			set popt="$popt $argv[$i1]"
			@ i1++
		end
		breaksw
	case -*:
		echo "Unknown option: $argv[$i]"
		exit 1
	default:
		set firstarg=$i
		breaksw
	endsw
	set i=$i1
end
set ecode=1
onintr quit
set td=`mktemp -d /tmp/xyz.XXXXXX`
if ( ! $firstarg ) then
	ra_xyze -r -u $popt > $td/stdin
	if ( $status ) goto quit
else
	set i=$firstarg
	while ( $i <= $#argv )
		ra_xyze -r -u $popt $argv[$i] $td/$argv[$i]:t
		if ( $status ) goto quit
		@ i++
	end
endif
ximage $xiargs $td/*
set ecode=$status
quit:
	rm -rf $td
	exit $ecode
notenough:
	echo "Missing arguments for $argv[$i] option"
	exit 1
