#!/bin/csh -f
# RCSid: $Id: phisto.csh,v 3.3 2003/02/22 02:07:27 greg Exp $
#
# Compute foveal histogram for picture set
#
set tf=/usr/tmp/ph$$
onintr quit
if ( $#argv == 0 ) then
	pfilt -1 -x 128 -y 128 -p 1 \
			| pvalue -o -h -H -d -b > $tf.dat
else
	rm -f $tf.dat
	foreach i ( $* )
		pfilt -1 -x 128 -y 128 -p 1 $i \
				| pvalue -o -h -H -d -b >> $tf.dat
		if ( $status ) exit 1
	end
endif
set Lmin=`total -l $tf.dat | rcalc -e 'L=$1*179;$1=if(L-1e-7,log10(L)-.01,-7)'`
set Lmax=`total -u $tf.dat | rcalc -e '$1=log10($1*179)+.01'`
rcalc -e 'L=$1*179;cond=L-1e-7;$1=log10(L)' $tf.dat \
	| histo $Lmin $Lmax 100
quit:
rm -f $tf.dat
