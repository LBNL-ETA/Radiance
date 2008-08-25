#!/bin/csh -f
# RCSid: $Id: phisto.csh,v 3.5 2008/08/25 04:50:32 greg Exp $
#
# Compute foveal histogram for picture set
#
set tf=`mktemp /tmp/phdat.XXXXX`
onintr quit
if ( $#argv == 0 ) then
	pfilt -1 -x 128 -y 128 -p 1 \
			| pvalue -o -h -H -d -b > $tf
else
	rm -f $tf
	foreach i ( $* )
		pfilt -1 -x 128 -y 128 -p 1 $i \
				| pvalue -o -h -H -d -b >> $tf
		if ( $status ) exit 1
	end
endif
set Lmin=`total -l $tf | rcalc -e 'L=$1*179;$1=if(L-1e-7,log10(L)-.01,-7)'`
set Lmax=`total -u $tf | rcalc -e '$1=log10($1*179)+.01'`
rcalc -e 'L=$1*179;cond=L-1e-7;$1=log10(L)' $tf \
	| histo $Lmin $Lmax 100
quit:
rm -f $tf
