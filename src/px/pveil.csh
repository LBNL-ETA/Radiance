#!/bin/csh -f
# RCSid: $Id: pveil.csh,v 3.3 2005/02/16 05:40:11 greg Exp $
#
# Add veiling glare to picture
#
if ($#argv != 1) then
	echo "Usage: $0 input.pic > output.pic"
	exit 1
endif
set ifile=$1
set td=/tmp
set gf=$td/av$$.gs
set cf=$td/av$$.cal
set tf=($gf $cf)
onintr quit
findglare -r 400 -c -p $ifile \
	| sed -e '1,/^BEGIN glare source$/d' -e '/^END glare source$/,$d' \
	> $gf
if ( -z $gf ) then
	cat $ifile
	goto quit
endif
( rcalc -e '$1=recno;$2=$1;$3=$2;$4=$3;$5=$4*$5' $gf \
	| tabfunc SDx SDy SDz I ; cat ) > $cf << '_EOF_'
N : I(0);
K : 9.2;	{ should be 9.6e-3/PI*(180/PI)^2 == 10.03 ? }
bound(a,x,b) : if(a-x, a, if(x-b, b, x));
Acos(x) : acos(bound(-1,x,1));
sq(x) : x*x;
mul(ct) : if(ct-cos(.5*PI/180), K/sq(.5), K/sq(180/PI)*ct/sq(Acos(ct)));
Dx1 = Dx(1); Dy1 = Dy(1); Dz1 = Dz(1);		{ minor optimization }
cosa(i) = SDx(i)*Dx1 + SDy(i)*Dy1 + SDz(i)*Dz1;
sum(i) = if(i-.5, mul(cosa(i))*I(i)+sum(i-1), 0);
veil = le(1)/WE * sum(N);
ro = ri(1) + veil;
go = gi(1) + veil;
bo = bi(1) + veil;
'_EOF_'
getinfo < $ifile | egrep '^((VIEW|EXPOSURE|PIXASPECT|PRIMARIES|COLORCORR)=|[^ ]*(rpict|rview|pinterp) )'
pcomb -f $cf $ifile
quit:
rm -f $tf
