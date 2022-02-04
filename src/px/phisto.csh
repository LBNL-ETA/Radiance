#!/bin/csh -f
# RCSid: $Id: phisto.csh,v 3.8 2022/02/04 20:11:49 greg Exp $
#
# Compute foveal histogram for picture set
#
set tf=`mktemp /tmp/phdat.XXXXX`
onintr quit
if ( $#argv == 0 ) then
	pfilt -1 -x 128 -y 128 -p 1 \
			| pvalue -O -h -H -df -b > $tf
else
	rm -f $tf
	foreach i ( $* )
		pfilt -1 -x 128 -y 128 -p 1 $i \
				| pvalue -O -h -H -df -b >> $tf
		if ( $status ) exit 1
	end
endif
set Lmin=`total -if -l $tf | rcalc -e 'L=$1*179;$1=if(L-1e-7,log10(L)-.01,-7)'`
set Lmax=`total -if -u $tf | rcalc -e '$1=log10($1*179)+.01'`
rcalc -if -e 'L=$1*179;cond=L-1e-7;$1=log10(L)' $tf \
	| histo $Lmin $Lmax 777
quit:
rm -f $tf
