#!/bin/csh -f
# RCSid $Id: pgblur.csh,v 1.1 2004/01/29 22:19:13 greg Exp $
#
# Apply Gaussian blur without resizing image
# More efficient than straight pfilt for large blurs
#
if ( $#argv != 3 ) then
	goto userr
endif
if ( "$1" != "-r" ) then
	goto userr
endif
if ( "$2" !~ [1-9]* ) then
	goto userr
endif
set rad = $2
set inp = "$3"
set reduc = `ev "floor($rad/1.8)"`
if ( $reduc <= 1 ) then
	exec pfilt -1 -r $rad $inp:q
endif
set filt = `ev "$rad/$reduc"`
set pr=`getinfo -d < $inp:q | sed 's/^-Y \([1-9][0-9]*\) +X \([1-9][0-9]*\)$/\2 \1/'`
pfilt -1 -x /$reduc -y /$reduc $inp:q \
	| pfilt -1 -r $filt -x $pr[1] -y $pr[2]
exit 0
userr:
echo Usage: "$0 -r radius input.pic"
exit 1
