#!/bin/csh -fe
# RCSid: $Id: falsecolor.csh,v 2.13 2005/11/12 06:14:51 greg Exp $
#
# Create false color image with legend
#
# Added user-definable legend 2004/01/20 Rob Guglielmetti

set td=/tmp/fc$$
onintr quit
set mult=179
set label=Nits
set scale=1000
set decades=0
# set redv='1.6*v-.6'
# set grnv='if(v-.375,1.6-1.6*v,8/3*v)'
# set bluv='1-8/3*v'
set redv='def_red(v)'
set grnv='def_grn(v)'
set bluv='def_blu(v)'
set ndivs=8
set picture='-'
set cpict=
set loff=0
set legwidth=100
set legheight=200
while ($#argv > 0)
	switch ($argv[1])
	case -lw:
		shift argv
		set legwidth="$argv[1]"
		breaksw
	case -lh:
		shift argv
		set legheight="$argv[1]"
		breaksw
	case -m:
		shift argv
		set mult="$argv[1]"
		breaksw
	case -s:
		shift argv
		set scale="$argv[1]"
		breaksw
	case -l:
		shift argv
		set label="$argv[1]"
		breaksw
	case -log:
		shift argv
		set decades=$argv[1]
		breaksw
	case -r:
		shift argv
		set redv="$argv[1]"
		breaksw
	case -g:
		shift argv
		set grnv="$argv[1]"
		breaksw
	case -b:
		shift argv
		set bluv="$argv[1]"
		breaksw
	case -i:
		shift argv
		set picture="$argv[1]"
		breaksw
	case -p:
		shift argv
		set cpict="$argv[1]"
		breaksw
	case -ip:
	case -pi:
		shift argv
		set picture="$argv[1]"
		set cpict="$argv[1]"
		breaksw
	case -cl:
		set docont=a
		set loff=12
		breaksw
	case -cb:
		set docont=b
		set loff=13
		breaksw
	case -e:
		set doextrem
		set needfile
		breaksw
	case -n:
		shift argv
		set ndivs="$argv[1]"
		breaksw
	default:
		echo bad option "'$argv[1]'"
		exit 1
	endsw
	shift argv
end
mkdir $td
cat > $td/pc0.cal <<_EOF_
PI : 3.14159265358979323846 ;
scale : $scale ;
mult : $mult ;
ndivs : $ndivs ;

or(a,b) : if(a,a,b);
EPS : 1e-7;
neq(a,b) : if(a-b-EPS,1,b-a-EPS);
btwn(a,x,b) : if(a-x,-1,b-x);
clip(x) : if(x-1,1,if(x,x,0));
frac(x) : x - floor(x);
boundary(a,b) : neq(floor(ndivs*a+.5),floor(ndivs*b+.5));

interp_arr2(i,x,f):(i+1-x)*f(i)+(x-i)*f(i+1);
interp_arr(x,f):if(x-1,if(f(0)-x,interp_arr2(floor(x),x,f),f(f(0))),f(1));
def_redp(i):select(i,0.18848,0.05468174,
0.00103547,0,7.449763e-06,0.0004390987,0.001367254,0.003076,
0.01376382,0.06085327,0.1712073,0.2889212,0.3404111,0.3803934,
0.4135148,0.4461151,0.4966744,0.6343959,0.7709373,0.9048474,1,
1,0.9863);
def_red(x):interp_arr(x/0.0454545+1,def_redp);
def_grnp(i):select(i,0.0009766,0,0.0008966244,
0.0264977,0.1256843,0.2865799,0.4251992,0.4775091,0.4406337,
0.3546771,0.2374256,0.1424751,0.09136063,0.0474073,0.0181899,
0.004181549,9.587855e-05,0,0.008718861,0.1025931,0.3106831,
0.6447838,0.9707);
def_grn(x):interp_arr(x/0.0454545+1,def_grnp);
def_blup(i):select(i,0.2666,0.3638662,0.4770437,
0.5131397,0.5363797,0.5193677,0.4091852,0.1732726,0.05331749,
0.0457587,0.06656252,0.06819438,0.04973373,0.02796644,
0.01131437,0.00280772,0.0002246617,0.0003480934,0.001044561,
0.001006845,0.002297005,0.01128842,0.02539);
def_blu(x):interp_arr(x/0.0454545+1,def_blup);

isconta = if(btwn(0,v,1),or(boundary(vleft,vright),boundary(vabove,vbelow)),-1);
iscontb = if(btwn(0,v,1),btwn(.4,frac(ndivs*v),.6),-1); 

ra = 0;
ga = 0;
ba = 0;

in = 1;

ro = if(in,clip($redv),ra);
go = if(in,clip($grnv),ga);
bo = if(in,clip($bluv),ba);
_EOF_
cat > $td/pc1.cal <<_EOF_
norm : mult/scale/le(1);

v = map(li(1)*norm);

vleft = map(li(1,-1,0)*norm);
vright = map(li(1,1,0)*norm);
vabove = map(li(1,0,1)*norm);
vbelow = map(li(1,0,-1)*norm);

map(x) = x;

ra = ri(nfiles);
ga = gi(nfiles);
ba = bi(nfiles);
_EOF_
set pc0args=(-f $td/pc0.cal)
set pc1args=(-f $td/pc1.cal)
if ($?docont) then
	set pc0args=($pc0args -e "in=iscont$docont")
endif
if ("$cpict" == "") then
	set pc1args=($pc1args -e 'ra=0;ga=0;ba=0')
else if ("$cpict" == "$picture") then
	set cpict=
endif
if ($?needfile && "$picture" == '-') then
	cat > $td/picture
	set picture=$td/picture
endif
if ("$decades" != "0") then
	set pc1args=($pc1args -e "map(x)=if(x-10^-$decades,log10(x)/$decades+1,0)")
	set imap="imap(y)=10^((y-1)*$decades)"
else
	set imap="imap(y)=y"
endif
if ( $legwidth > 20 && $legheight > 40 ) then
pcomb $pc0args -e 'v=(y+.5)/yres;vleft=v;vright=v' \
		-e 'vbelow=(y-.5)/yres;vabove=(y+1.5)/yres' \
		-x $legwidth -y $legheight > $td/scol.pic
( echo "$label"; cnt $ndivs \
		| rcalc -e '$1='"($scale)*imap(($ndivs-.5-"'$1'")/$ndivs)" \
		-e "$imap" | sed -e 's/\(\.[0-9][0-9][0-9]\)[0-9]*/\1/' ) \
	| psign -s -.15 -cf 1 1 1 -cb 0 0 0 \
		-h `ev "floor($legheight/$ndivs+.5)"` > $td/slab.pic
else
	set legwidth=0
	set legheight=0
	(echo "" ; echo "-Y 1 +X 1" ; echo "aaa" ) > $td/scol.pic
	cp $td/scol.pic $td/slab.pic
endif
if ( $?doextrem ) then
	pextrem -o $picture > $td/extrema
	set minpos=`sed 2d $td/extrema | rcalc -e '$2=$2;$1=$1+'"$legwidth"`
	set minval=`rcalc -e '$1=($3*.27+$4*.67+$5*.06)*'"$mult" $td/extrema | sed -e 2d -e 's/\(\.[0-9][0-9][0-9]\)[0-9]*/\1/'`
	set maxpos=`sed 1d $td/extrema | rcalc -e '$2=$2;$1=$1+'"$legwidth"`
	set maxval=`rcalc -e '$1=($3*.27+$4*.67+$5*.06)*'"$mult" $td/extrema | sed -e 1d -e 's/\(\.[0-9][0-9][0-9]\)[0-9]*/\1/'`
	psign -s -.15 -a 2 -h 16 $minval > $td/minv.pic
	psign -s -.15 -a 2 -h 16 $maxval > $td/maxv.pic
	pcomb $pc0args $pc1args $picture $cpict \
		| pcompos $td/scol.pic 0 0 \
			+t .1 "\!pcomb -e 'lo=1-gi(1)' $td/slab.pic" \
			`ev 2 $loff-1` -t .5 $td/slab.pic 0 $loff \
		  - $legwidth 0 $td/minv.pic $minpos $td/maxv.pic $maxpos
else
	pcomb $pc0args $pc1args $picture $cpict \
		| pcompos $td/scol.pic 0 0 \
			+t .1 "\!pcomb -e 'lo=1-gi(1)' $td/slab.pic" \
			`ev 2 $loff-1` -t .5 $td/slab.pic 0 $loff \
			- $legwidth 0
endif
quit:
rm -rf $td
