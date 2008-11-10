#!/bin/csh -f
# RCSid: $Id: psquish.csh,v 3.5 2008/11/10 19:08:19 greg Exp $
set Lmin=.0001		# minimum visible world luminance
set Ldmin=1		# minimum display luminance
set Ldmax=100		# maximum display luminance
set nsteps=100		# number of steps in perceptual histogram
set cvratio=0.05	# fraction of pixels to ignore in envelope clipping
set td=/tmp
set tf0=$td/tf$$
set tf1=$td/hist$$
set tf1b=$td/hist$$.diff
set tf2=$td/cumt$$
set tf3=$td/histeq$$.cal
set tf4=$td/cf$$.cal
set tf=($tf0 $tf1 $tf1b $tf2 $tf3 $tf4)
if ( "$argv[1]" == "-a" ) then
	set adaptive
	shift argv
endif
if ( $#argv != 1 ) then
	echo "Usage: $0 [-a] input.hdr > output.hdr"
	exit 1
endif
set ifile=$1
set ibase=$ifile:t
if ( "$ibase" =~ *.hdr ) set ibase=$ibase:r
set ibase=$ibase:t
onintr quit
pextrem -o $ifile > $tf0
set Lmin=`rcalc -e 'L=179*(.265*$3+.67*$4+.065*$5)' -e 'cond=1.5-recno;$1=if('L-$Lmin,L,$Lmin')' $tf0`
set Lmax=`rcalc -e 'L=179*(.265*$3+.67*$4+.065*$5)' -e 'cond=recno-1.5;$1=if('$Ldmax-L,$Ldmax,L')' $tf0`
cat > $tf3 << _EOF_
min(a,b) : if(a-b, b, a);
WE : 179;			{ Radiance white luminous efficacy }
Lmin : $Lmin ;			{ minimum visible luminance }
Lmax : $Lmax ;			{ maximum picture luminance }
Ldmin : $Ldmin ;		{ minimum output luminance }
Ldmax : $Ldmax ;		{ maximum output luminance }
Stepsiz : (Bl(Lmax)-Bl(Lmin))/ $nsteps ;	{ brightness step size }
			{ Logarithmic brightness function }
Bl(L) : log(L);
Lb(B) : exp(B);
BLw(Lw) : Bl(Ldmin) + (Bl(Ldmax)-Bl(Ldmin))*cf(Bl(Lw));
			{ first derivative functions }
Bl1(L) : 1/L;
Lb1(B) : exp(B);
			{ histogram equalization function }
lin = li(1);
Lw = WE/le(1) * lin;
Lout = Lb(BLw(Lw));
mult = if(Lw-Lmin, (Lout-Ldmin)/(Ldmax-Ldmin)/lin, 0) ;
_EOF_
if ( $?adaptive ) then
	cat >> $tf3 << _EOF_
			{ Ferwerda contrast sensitivity function }
				{ log10 of cone threshold luminance }
ltp(lLa) : if(-2.6 - lLa, -.72, if(lLa - 1.9, lLa - 1.255,
		(.249*lLa + .65)^2.7 - .72));
				{ log10 of rod threshold luminance }
lts(lLa) : if(-3.94 - lLa, -2.86, if(lLa - -1.44, lLa - .395,
		(.405*lLa + 1.6)^2.18 - 2.86));
				{ threshold is minimum of rods and cones }
ldL2(lLa) : min(ltp(lLa),lts(lLa));
dL(La) : 10^ldL2(log10(La));
			{ derivative clamping function }
clamp2(L, bLw) : dL(Lb(bLw))/dL(L)/Lb1(bLw)/(Bl(Ldmax)-Bl(Ldmin))/Bl1(L);
clamp(L) : clamp2(L, BLw(L));
			{ shift direction for histogram }
shiftdir(B) : if(cf(B) - (B - Bl(Lmin))/(Bl(Ldmax) - Bl(Lmin)), 1, -1);
			{ Scotopic/Photopic color adjustment }
sL(r,g,b) : .062*r + .608*g + .330*b;	{ approx. scotopic brightness }
BotMesopic : 10^-2.25;		{ top of scotopic range }
TopMesopic : 10^0.75;		{ bottom of photopic range }
incolor = if(Lw-TopMesopic, 1, if(BotMesopic-Lw, 0,
		(Lw-BotMesopic)/(TopMesopic-BotMesopic)));
slf = (1 - incolor)*sL(ri(1),gi(1),bi(1));
ro = mult*(incolor*ri(1) + slf);
go = mult*(incolor*gi(1) + slf);
bo = mult*(incolor*bi(1) + slf);
_EOF_
else
	cat >> $tf3 << _EOF_
			{ derivative clamping function }
clamp2(L, bLw) : Lb(bLw)/L/Lb1(bLw)/(Bl(Ldmax)-Bl(Ldmin))/Bl1(L);
clamp(L) : clamp2(L, BLw(L));
			{ shift direction for histogram }
shiftdir(B) : -1;
ro = mult*ri(1);
go = mult*gi(1);
bo = mult*bi(1);
_EOF_
endif
# Compute brightness histogram
pfilt -1 -p 1 -x 128 -y 128 $ifile | pvalue -o -b -d -h -H \
	| rcalc -f $tf3 -e 'Lw=WE*$1;$1=if(Lw-Lmin,Bl(Lw),-1)' \
	| histo `ev "log($Lmin)" "log($Lmax)"` $nsteps > $tf1
# Clamp frequency distribution
set totcount=`sed 's/^.*[ 	]//' $tf1 | total`
set margin=`ev "floor($totcount*$cvratio+.5)"`
while ( 1 )
	# Compute mapping function
	sed 's/^.*[ 	]//' $tf1 | total -1 -r \
		| rcalc -e '$1=$1/'$totcount | rlam $tf1 - \
		| tabfunc -i 0 cf > $tf4
	# Compute difference with visible envelope
	rcalc -f $tf4 -f $tf3 -e "T:$totcount*Stepsiz" \
			-e 'clfq=floor(T*clamp(Lb($1))+.5)' \
			-e '$1=$2-clfq;$2=shiftdir($1)' $tf1 > $tf1b
	if (`sed 's/[	 ].*$//' $tf1b | total` >= 0) then
		# Nothing visible? -- just normalize
		pfilt -1 -e `pextrem $ifile | rcalc -e 'cond=recno-1.5;$1=1/(.265*$3+.67*$4+.065*$5)'` $ifile
		goto quit
	endif
	# Check to see if we're close enough
	if (`rcalc -e '$1=if($1,$1,0)' $tf1b | total` <= $margin) break
	# Squash frequency distribution
	set diffs=(`sed 's/[ 	].*$//' $tf1b`)
	set shftd=(`sed 's/^.*[ 	]//' $tf1b`)
	while ( 1 )
		set maxi=0
		set maxd=0
		set i=$nsteps
		while ( $i > 0 )
			if ( $diffs[$i] > $maxd ) then
				set maxd=$diffs[$i]
				set maxi=$i
			endif
			@ i--
		end
		if ( $maxd == 0 ) break
		set i=0
	tryagain:
		set r=$shftd[$maxi]
		while ( $i == 0 )
			@ t= $maxi + $r
			if ( $t < 1 || $t > $nsteps ) then
				@ shftd[$maxi]= -($shftd[$maxi])
				goto tryagain
			endif
			if ( $diffs[$t] < 0 ) set i=$t
			@ r+= $shftd[$maxi]
		end
		if ( $diffs[$i] <= -$diffs[$maxi] ) then
			@ diffs[$i]+= $diffs[$maxi]
			set diffs[$maxi]=0
		else
			@ diffs[$maxi]+= $diffs[$i]
			set diffs[$i]=0
		endif
	end
	# Mung histogram
	echo $diffs | tr ' ' '\012' | rlam $tf1 - \
		| rcalc -f $tf4 -f $tf3 -e "T:$totcount*Stepsiz" \
			-e 'clfq=floor(T*clamp(Lb($1))+.5)' \
			-e '$1=$1;$2=$3+clfq' > $tf1b
	mv -f $tf1b $tf1
end
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
		-e 'lx=$1/99*(log10(Lmax)-log10(Lmin))+log10(Lmin)' \
		>> ${ibase}_brmap.plt
	if ( $?DISPLAY ) then
		bgraph ${ibase}_histo.plt ${ibase}_brmap.plt | x11meta &
	endif
endif
# Map our picture
getinfo < $ifile | egrep '^((VIEW|PIXASPECT|PRIMARIES)=|[^ ]*(rpict|rview|pinterp) )'
pcomb -f $tf4 -f $tf3 $ifile
quit:
rm -f $tf
