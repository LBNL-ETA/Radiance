#!/bin/csh -f
# RCSid: $Id: illumcal.csh,v 1.1 2003/02/22 02:07:21 greg Exp $
#
# Compute color characteristics of light sources
#
# Each input file should contain evenly-spaced pairs of wavelength (nm) and
# power values (watts/nm), one per line.
#
if ( $#argv < 1 ) then
	echo Usage: $0 illum.dat ..
	exit 1
endif
set cal = .
foreach illum ($*)
set spc=(`sed -e 1d -e 's/^[ 	]*\([1-9][0-9]*\)[ 	].*$/\1/' -e 3q $illum`)
rcalc -f $cal/cieresp.cal -f $cal/stdrefl.cal -e "intvl=abs($spc[2]-$spc[1])" \
	-f $cal/conv1.cal $illum | total >> /tmp/il$$.dat
end
rcalc -f $cal/conv2.cal -f $cal/cct.cal -f $cal/cri.cal \
	-o $cal/illum.fmt /tmp/il$$.dat
rm -f /tmp/il$$.dat
