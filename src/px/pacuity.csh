#!/bin/csh -f
# RCSid: $Id: pacuity.csh,v 3.6 2008/11/10 19:08:19 greg Exp $
#
# Adjust picture acuity according to human visual abilities
#
if ($#argv != 1) then
	echo "Usage: $0 input.hdr > output.hdr"
	exit 1
endif
set td=/tmp
set tfc1=$td/ac$$.cal
set tf=($td/c{1,2,4,8,16,32}d$$.hdr $tfc1)
set ifile=$1
onintr quit
tabfunc -i acuity > $tfc1 << EOF
# Log10 luminance vs. visual acuity in cycles/degree
-2.804	2.09
-2.363	3.28
-2.076	3.79
-1.792	4.39
-1.343	6.11
-1.057	8.83
-0.773	10.94
-0.371	18.66
-0.084	23.88
0.2	31.05
0.595	37.42
0.882	37.68
1.166	41.60
1.558	43.16
1.845	45.30
2.129	47.00
2.577	48.43
2.864	48.32
3.148	51.06
3.550	51.09
EOF
set pres=(`getinfo -d < $ifile | sed 's/^-Y \([1-9][0-9]*\) +X \([1-9][0-9]*\)$/\2 \1/'`)
set vp=`vwright V < $ifile`
set aext=(`pextrem -o $ifile | rcalc -f $tfc1 -e 'max(a,b):if(a-b,a,b);$1=acuity(log10(max(179*(.265*$3+.67*$4+.065*$5),1e-4)))'`)
( rcalc -e "$vp" -e "A:3438*sqrt(Vhn/$pres[1]*Vvn/$pres[2])" \
	-e 'f=60/A/2/$1;cond=if(1.5-$1,1,if(1-f,-1,if($1-'"$aext[2]"',-1,$1-'"$aext[1])))" \
	-o 'pfilt -1 -r 2 -x /${f} -y /${f} '"$ifile | pfilt -1 -r 1 -x $pres[1] -y $pres[2] > $td/"'c${$1}d'$$.hdr \
	| csh -f ) << EOF
1
2
4
8
16
32
EOF
cat >> $tfc1 << _EOF_
max(a,b) : if(a-b, a, b);
target_acuity = acuity(log10(max(WE/le(1)*li(1),1e-4)));
findfuzzy(i) = if(target_acuity-picture_acuity(i),i,if(i-1.5,findfuzzy(i-1),1));
fuzzy_picture = findfuzzy(nfiles-1);
clear_picture = fuzzy_picture + 1;
clarity_ex = (target_acuity-picture_acuity(fuzzy_picture)) /
		(picture_acuity(clear_picture)-picture_acuity(fuzzy_picture));
clarity = if(clarity_ex-1, 1, if(-clarity_ex, 0, clarity_ex));
ro = clarity*ri(clear_picture) + (1-clarity)*ri(fuzzy_picture);
go = clarity*gi(clear_picture) + (1-clarity)*gi(fuzzy_picture);
bo = clarity*bi(clear_picture) + (1-clarity)*bi(fuzzy_picture);
picture_acuity(n) : select(n,1,
_EOF_
set pf=($td/c1d$$.hdr)
foreach i (2 4 8 16 32)
	if ( -f $td/c${i}d$$.hdr ) then
		set pf=( $pf $td/c${i}d$$.hdr )
		echo -n "$i," >> $tfc1
	endif
end
set pf=( $pf $ifile )
rcalc -n -e "$vp" -e "A:3879*sqrt(Vhn/$pres[1]*Vvn/$pres[2])" \
	-o '${60/A});' >> $tfc1
getinfo < $ifile | egrep '^((VIEW|EXPOSURE|PIXASPECT|PRIMARIES|COLORCORR)=|[^ ]*(rpict|rview|pinterp) )'
pcomb -f $tfc1 $pf
quit:
rm -f $tf
