#!/bin/csh -f
# RCSid: $Id: vlpic.csh,v 3.2 2005/01/18 03:59:41 greg Exp $
#
# Compute falsecolor image of visibility level
# using the wacky formulas of Werner Adrian.
#
set age=40			# default age (years)
set tim=.2			# default time (seconds)
while ($#argv > 1)
	switch ($argv[1])
	case -a:
		shift argv
		set age="$argv[1]"
		breaksw
	case -t:
		shift argv
		set tim="$argv[1]"
		breaksw
	default:
		echo bad option "'$argv[1]'"
		exit 1
	endsw
	shift argv
end
set tc=/usr/tmp/vl$$.cal
set tp1=/usr/tmp/vl$$r1.pic
set tp2=/usr/tmp/vl$$r2.pic
set tp4=/usr/tmp/vl$$r4.pic
set tp8=/usr/tmp/vl$$r8.pic
set tf=($tc $tp1 $tp2 $tp4 $tp8)
set inpic=$argv[1]
onintr quit
set pr=(`getinfo -d < $inpic | sed 's/^-Y \([1-9][0-9]*\) +X \([1-9][0-9]*\)$/\2 \1/'`)
# ( vwright V < $inpic ; cat ) > $tc << _EOF_
cat > $tc << _EOF_
{ A : 3438 * sqrt(Vhn/xmax*Vvn/ymax);	{ pixel size (in minutes) } }
PI : 3.14159265358979323846;
A = 180*60/PI * sqrt(S(1)/PI);
sq(x) : x*x;
max(a,b) : if(a-b, a, b);
Lv : 0;		{ veiling luminance {temporarily zero} }
		{ find direction of maximum contrast }
xoff(d) = select(d, 1, 1, 0, -1, -1, -1, 0, 1);
yoff(d) = select(d, 0, 1, 1, 1, 0, -1, -1, -1);
cf(lf, lb) = (lf - lb)/(lb + Lv);
contrast(d) = cf(li(1,xoff(d),yoff(d)), li(1,-xoff(d),-yoff(d)));
cwin(d1, d2) = if(contrast(d1) - contrast(d2), d1, d2);
bestdir(d) = if(d-1.5, cwin(d, bestdir(d-1)), d);
maxdir = bestdir(8);
Lt = WE*li(1,xoff(maxdir),yoff(maxdir)) + Lv;
La = WE*li(1,-xoff(maxdir),-yoff(maxdir)) + Lv;
LLa = log10(La);
F = if(La-.6, log10(4.1925) + LLa*.1556 + .1684*La^.5867,
	if(La-.00418, 10^(-.072 + LLa*(.3372 + LLa*.0866)),
		10^(.028 + .173*LLa)));
L = if(La-.6, .05946*La^.466,
	if(La-.00418, 10^(-1.256 + .319*LLa),
		10^(-.891 + LLa*(.5275 + LLa*.0227))));
LAt = log10(A) + .523;
AA = .36 - .0972*LAt*LAt/(LAt*(LAt - 2.513) + 2.7895);
LLat = LLa + 6;
AL = .355 - .1217*LLat*LLat/(LLat*(LLat - 10.4) + 52.28);
AZ = sqrt(AA*AA + AL*AL)/2.1;
DL1 = 2.6*sq(F/A + L);
M = if(La-.004, 10^-(10^-(if(La-.1,.125,.075)*sq(LLa+1) + .0245)), 0);
TGB = .6*La^-.1488;
FCP = 1 - M*A^-TGB/(2.4*(DL1*(AZ+2)/2));
DL2 = DL1*(AZ + T)/T;
FA = if(Age-64, sq(Age-56.6)/116.3 + 1.43, sq(Age-19)/2160 + .99);
DL3 = DL2 * FA;
DL4 = if(La-Lt, DL3*FCP, DL3);
lo = (Lt - La)/DL4;		{ Output VL }
_EOF_
pcomb -w -e "Age:$age;T:$tim" -f $tc -o $inpic > $tp1
pfilt -1 -x /2 -y /2 $inpic \
	| pcomb -w -e "Age:$age;T:$tim" -f $tc -o - \
	| pfilt -1 -r 1 -x $pr[1] -y $pr[2] > $tp2
pfilt -1 -x /4 -y /4 $inpic \
	| pcomb -w -e "Age:$age;T:$tim" -f $tc -o - \
	| pfilt -1 -r 1 -x $pr[1] -y $pr[2] > $tp4
pfilt -1 -x /8 -y /8 $inpic \
	| pcomb -w -e "Age:$age;T:$tim" -f $tc -o - \
	| pfilt -1 -r 1 -x $pr[1] -y $pr[2] > $tp8
pcomb -e 'max(a,b):if(a-b,a,b)' -e 'lo=max(gi(1),max(gi(2),max(gi(3),gi(4))))' \
	$tp1 $tp2 $tp4 $tp8 \
	| falsecolor -l VL -s 56 -m 1 -r 'if(v-1/8,1.6*8/7*(v-1/8)-.6,0)' \
		-g 'if(v-1/8,if(v-.453125,1.6-1.6*8/7*(v-1/8),8/3*8/7*(v-1/8)),0)' \
		-b 'if(v-1/8,1-8/3*8/7*(v-1/8),0)'
quit:
rm -f $tf
