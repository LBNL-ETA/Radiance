#!/bin/csh -f
# RCSid: $Id: phisteq.csh,v 3.2 2004/01/01 19:44:07 greg Exp $
set Ldmin=1		# minimum display luminance
set Ldmax=100		# maximum display luminance
set nsteps=100		# number of steps in perceptual histogram
set cvratio=0.08	# fraction of pixels to ignore in envelope clipping
set td=/usr/tmp
set tf1=$td/hist$$
set tf1b=$td/hist$$.new
set tf2=$td/cumt$$
set tf3=$td/histeq$$.cal
set tf4=$td/cf$$.cal
set tf=($tf1 $tf1b $tf2 $tf3 $tf4)
if ( $#argv != 1 ) then
	echo "Usage: $0 input.pic > output.pic"
	exit 1
endif
set ifile=$1
set ibase=$ifile:t
if ( "$ibase" =~ *.pic ) set ibase=$ibase:r
set ibase=$ibase:t
onintr quit
cat > $tf3 << _EOF_
WE : 179;			{ Radiance white luminous efficacy }
Lmin : .0001;			{ minimum allowed luminance }
Ldmin : $Ldmin ;		{ minimum output luminance }
Ldmax : $Ldmax ;		{ maximum output luminance }
Stepsiz : 1/ $nsteps ;		{ brightness step size }
			{ Daly local amplitude nonlinearity formulae }
sq(x) : x*x;
c1 : 12.6;
b : .63;
Bl(L) : L / (L + (c1*L)^b);
Lb(B) : (c1^b*B/(1-B))^(1/(1-b));
BLw(Lw) : Bl(Ldmin) + (Bl(Ldmax)-Bl(Ldmin))*cf(Bl(Lw));
			{ first derivative functions }
Bl1(L) : (c1*L)^b*(1-b)/sq(L + (c1*L)^b);
Lb1(B) : c1^b/(1-b)/sq(1-B) * (c1^b*B/(1-B))^(b/(1-b));
			{ derivative clamping function }
clamp2(L, aLw) : Lb(aLw) / L / Lb1(aLw) / (Bl(Ldmax)-Bl(Ldmin)) / Bl1(L);
clamp(L) : clamp2(L, BLw(L));
			{ histogram equalization function }
lin = li(1);
Lw = WE/le(1) * lin;
Lout = Lb(BLw(Lw));
mult = if(Lw-Lmin, (Lout-Ldmin)/(Ldmax-Ldmin)/lin, 0) ;
ro = mult * ri(1);
go = mult * gi(1);
bo = mult * bi(1);
_EOF_
# Compute brightness histogram
pfilt -1 -p 1 -x 128 -y 128 $ifile | pvalue -o -b -d -h -H \
	| rcalc -f $tf3 -e 'Lw=WE*$1;$1=if(Lw-Lmin,Bl(Lw),-1)' \
	| histo 0 1 $nsteps | sed '/[ 	]0$/d' > $tf1
# Clamp frequency distribution
set totcount=`sed 's/^.*[ 	]//' $tf1 | total`
set tst=1
while ( $totcount > 0 )
	sed 's/^.*[ 	]//' $tf1 | total -1 -r \
		| rcalc -e '$1=$1/'$totcount | rlam $tf1 - \
		| tabfunc -i 0 cf > $tf4
	if ( $tst <= 0 ) break
	rcalc -f $tf4 -f $tf3 -e "T:$totcount*Stepsiz" \
			-e 'clfq=floor(T*clamp(Lb($1)))' \
			-e '$1=$1;$2=if($2-clfq,clfq,$2)' $tf1 > $tf1b
	set newtot=`sed 's/^.*[ 	]//' $tf1b | total`
	set tst=`ev "floor((1-$cvratio)*$totcount)-$newtot"`
	mv -f $tf1b $tf1
	set totcount=$newtot
end
if ( $totcount < 1 ) then
	# Fits in display range nicely already -- just normalize
	pfilt -1 -e `pextrem $ifile | rcalc -e 'cond=recno-1.5;$1=1/(.265*$3+.67*$4+.065*$5)'` $ifile
else
	# Plot the mapping function if we are in debug mode
	if ( $?DEBUG ) then
		cat > ${ibase}_histo.plt << _EOF_
include=curve.plt
title="Brightness Frequency Distribution"
subtitle= $ibase
ymin=0
xlabel="Perceptual Brightness B(Lw)"
ylabel="Frequency Count"
Alabel="Histogram"
Alintype=0
Blabel="Envelope"
Bsymsize=0
Adata=
_EOF_
		(cat $tf1; echo \;; echo Bdata=) >> ${ibase}_histo.plt
		rcalc -f $tf4 -f $tf3 -e "T:$totcount*Stepsiz" \
			-e '$1=$1;$2=T*clamp(Lb($1))' $tf1 \
			>> ${ibase}_histo.plt
		cat > ${ibase}_brmap.plt << _EOF_
include=line.plt
title="Brightness Mapping Function"
subtitle= $ibase
xlabel="World Luminance (log cd/m^2)"
ylabel="Display Luminance (cd/m^2)"
ymax= $Ldmax
Adata=
_EOF_
		cnt 100 | rcalc -f $tf4 -f $tf3 -e '$1=lx;$2=Lb(BLw(10^lx))' \
			-e Lmin:Lb\(`sed -n '1s/[ 	].*$//p' $tf1`\) \
			-e Lmax:Lb\(`sed -n '$s/[ 	].*$//p' $tf1`\) \
			-e 'lx=$1/99*(log10(Lmax)-log10(Lmin))+log10(Lmin)' \
			>> ${ibase}_brmap.plt
		if ( $?DISPLAY ) then
			bgraph ${ibase}_histo.plt ${ibase}_brmap.plt | x11meta &
		endif
	endif
	# Map our picture
	pcomb -f $tf4 -f $tf3 $ifile
endif
quit:
rm -f $tf
