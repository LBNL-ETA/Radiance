#!/bin/csh -f
# RCSid $Id: pbilat.csh,v 1.2 2008/08/25 04:50:32 greg Exp $
#
# Bilateral Filter (fixed parameters for now)
#
set bfrac=0.02
set sigma=0.4
set vmin=1e-7
if ( $#argv != 1 ) then
	goto userr
endif
set inp="$1"
set rad=`getinfo -d < $inp:q | rcalc -i '-Y ${yr} +X ${xr}' -e '$1=sqrt(xr*yr)*'$bfrac`
set extrem=`pextrem -o $inp:q | rcalc -e "vmin:$vmin" -e 'max(a,b):if(a-b,a,b);$1=max($3*.265+$4*.670+$5*.065,vmin)'`
set nseg=`ev "ceil(log10($extrem[2]/$extrem[1])/$sigma)"`
if ( $nseg > 20 ) set nseg=20
onintr done
set tdir=`mktemp -d /tmp/blf.XXXXXX`
set gfunc="sq(x):x*x;gfunc(x):if(sq(x)-sq($sigma),0,sq(1-sq(x/$sigma)))"
set i=0
set imglist=()
while ( $i <= $nseg )
	set intens=`ev "$extrem[1]*($extrem[2]/$extrem[1])^($i/$nseg)"`
	pcomb -e $gfunc:q -e "vmin:$vmin" \
			-e 'max(a,b):if(a-b,a,b);l1=max(li(1),vmin)' \
			-e "lo=gfunc(log10(l1)-log10($intens))" \
			-o $inp > $tdir/gimg.pic
	pcomb -e 'sf=gi(2);ro=sf*ri(1);go=sf*gi(1);bo=sf*bi(1)' \
			-o $inp $tdir/gimg.pic > $tdir/g_p.pic
	pgblur -r $rad $tdir/g_p.pic \
		| pcomb -e 'sf=if(gi(2)-1e-6,1/gi(2),1e6)' \
			-e 'ro=sf*ri(1);go=sf*gi(1);bo=sf*bi(1)' \
			- "\!pgblur -r $rad $tdir/gimg.pic" \
		> $tdir/gimg$i.pic
	set imglist=($imglist $tdir/gimg$i.pic)
	@ i++
end
cat > $tdir/interp.cal << _EOF_
{ Segmented intensity image interpolation }
max(a,b) : if(a-b, a, b);
NSEG : $nseg;
l1 = max(li(1), $vmin);
val = NSEG / log($extrem[2]/$extrem[1]) * (log(l1) - log($extrem[1]));
iv = floor(val);
ival = if(iv - (NSEG-1), NSEG-1, max(iv, 0));
xval = val - ival;
ro = (1-xval)*ri(ival+2) + xval*ri(ival+3);
go = (1-xval)*gi(ival+2) + xval*gi(ival+3);
bo = (1-xval)*bi(ival+2) + xval*bi(ival+3);
_EOF_
pcomb -h -f $tdir/interp.cal -o $inp:q $imglist
done:
rm -rf $tdir
exit 0
userr:
echo Usage: $0 input.pic
