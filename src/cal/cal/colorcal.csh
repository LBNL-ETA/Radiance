#!/bin/csh -f
# RCSid: $Id: colorcal.csh,v 1.1 2003/02/22 02:07:21 greg Exp $
#
# Compute CIE chromaticities from spectral reflectance data
#
if ( $#argv < 1 ) goto userr
set cal = .
if ( $argv[1] == "-i" ) then
	if ( $#argv < 3 ) goto userr
	shift argv
	set illum=$argv[1]
	shift argv
	foreach r ( $argv[*] )
		tabfunc -i rf < $r > /tmp/rf$$.cal
		rcalc -f $cal/cieresp.cal -f /tmp/rf$$.cal \
			-e 'r=rf($1);ty=$2*triy($1)' \
			-e '$1=ty;$2=$2*r*trix($1);$3=r*ty' \
			-e '$4=$2*r*triz($1)' \
			-e 'cond=if($1-359,831-$1,-1)' \
			$illum | total -m >> /tmp/rc$$.dat
	end
	rm -f /tmp/rf$$.cal
else
	foreach r ( $argv[*] )
		rcalc -f $cal/cieresp.cal -e 'ty=triy($1);$1=ty' \
		-e '$2=$2*trix($1);$3=$2*ty;$4=$2*triz($1)' \
		-e 'cond=if($1-359,831-$1,-1)' $r \
			| total -m >> /tmp/rc$$.dat
	end
endif
rcalc -e 'X=$2/$1;Y=$3/$1;Z=$4/$1' \
	-e 'x=X/(X+Y+Z);y=Y/(X+Y+Z);u=4*X/(X+15*Y+3*Z);v=9*Y/(X+15*Y+3*Z)' \
	-o $cal/color.fmt /tmp/rc$$.dat
rm -f /tmp/rc$$.dat
exit 0
userr:
echo "Usage: $0 [-i illum.dat] refl.dat .."
exit 1
