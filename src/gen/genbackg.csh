#!/bin/csh -f
# RCSid: $Id: genbackg.csh,v 2.1 2003/02/22 02:07:23 greg Exp $
#
# Generate height field to surround a mesh with an irregular border
#
#	4/16/2002	Greg Ward
#
if ( $#argv < 3 ) goto userr
# Get arguments
set mat = $1
set nam = $2
set mesh = $3
# Get options
set smooth=""
set i = 4
while ( $i <= $#argv )
	switch ($argv[$i])
	case -c:
		@ i++ ; set a=$i
		@ i++ ; set b=$i
		@ i++ ; set c=$i
		set c00=($argv[$a] $argv[$b] $argv[$c])
		@ i++ ; set a=$i
		@ i++ ; set b=$i
		@ i++ ; set c=$i
		set c01=($argv[$a] $argv[$b] $argv[$c])
		@ i++ ; set a=$i
		@ i++ ; set b=$i
		@ i++ ; set c=$i
		set c11=($argv[$a] $argv[$b] $argv[$c])
		@ i++ ; set a=$i
		@ i++ ; set b=$i
		@ i++ ; set c=$i
		set c10=($argv[$a] $argv[$b] $argv[$c])
		unset a b c
		breaksw
	case -r:
		@ i++
		set step = $argv[$i]
		breaksw
	case -b:
		@ i++
		set bord = $argv[$i]
		breaksw
	case -s:
		set smooth="-s"
		breaksw
	default:
		goto userr
	endsw
	@ i++
end
# Get corner points
set lim=(`getbbox -h $mesh`)
if (! $?bord ) then
	set bord=`ev "($lim[2]-$lim[1]+$lim[4]-$lim[3])/20"`
endif
if (! $?c00 ) then
	set c00=(`ev "$lim[1]-$bord" "$lim[3]-$bord" "($lim[5]+$lim[6])/2"`)
	set c01=(`ev "$lim[2]+$bord" "$lim[3]-$bord" "($lim[5]+$lim[6])/2"`)
	set c10=(`ev "$lim[1]-$bord" "$lim[4]+$bord" "($lim[5]+$lim[6])/2"`)
	set c11=(`ev "$lim[2]+$bord" "$lim[4]+$bord" "($lim[5]+$lim[6])/2"`)
endif
# Get mesh resolution
if (! $?step ) then
	set step=`ev "sqrt(($c11[1]-$c00[1])*($c11[1]-$c00[1])+($c11[2]-$c00[2])*($c11[2]-$c00[2]))/70"`
endif
set res=(`ev "ceil(sqrt(($c01[1]-$c00[1])*($c01[1]-$c00[1])+($c01[2]-$c00[2])*($c01[2]-$c00[2]))/$step)" "ceil(sqrt(($c10[1]-$c00[1])*($c10[1]-$c00[1])+($c10[2]-$c00[2])*($c10[2]-$c00[2]))/$step)"`)
cat > /tmp/edge$$.rad << _EOF_
void sphere c00
0 0 4
$c00 .001
void sphere c01
0 0 4
$c01 .001
void sphere c10
0 0 4
$c10 .001
void sphere c11
0 0 4
$c11 .001
void cylinder ex0
0 0 7
$c00 $c10 .001
void cylinder ex1
0 0 7
$c01 $c11 .001
void cylinder ey0
0 0 7
$c00 $c01 .001
void cylinder ey1
0 0 7
$c10 $c11 .001
_EOF_
cat > /tmp/mesh$$.cal << _EOF_
lerp(x,y0,y1):(1-x)*y0+x*y1;
and(a,b):if(a,b,-1);
or(a,b):if(a,1,b);
not(a):if(a,-1,1);
max(a,b):if(a-b,a,b);
min(a,b):if(a-b,b,a);
EPS:1e-7;
XR:min(32,floor($res[1]/2));
YR:8;
eq(a,b):and(a-b+EPS,b-a+EPS);
sumfun2x(f,y,x0,x1):if(x1-x0+EPS,f(x0,y)+sumfun2x(f,y,x0+1,x1),0);
sumfun2(f,x0,x1,y0,y1):if(y1-y0+EPS,
	sumfun2x(f,y0,x0,x1)+sumfun2(f,x0,x1,y0+1,y1), 0);
xpos(xf,yf)=lerp(xf,lerp(yf,$c00[1],$c10[1]),lerp(yf,$c01[1],$c11[1]));
ypos(xf,yf)=lerp(xf,lerp(yf,$c00[2],$c10[2]),lerp(yf,$c01[2],$c11[2]));
onedge=or(or(eq(x,0),eq(x,xmax-1)),or(eq(y,0),eq(y,ymax-1)));
z0:$lim[5]-$bord;
height(xo,yo)=z0+gi(1,xo,yo);
maxwt(xo,yo):if(and(eq(xo,0),eq(yo,0)),100,1/(xo*xo+yo*yo));
wt(xo,yo)=if(height(xo,yo)-1e9, 0, maxwt(xo,yo));
wtv(xo,yo)=wt(xo,yo)*height(xo,yo);
wts=sumfun2(wt,-XR,XR,-YR,YR);
sum=sumfun2(wtv,-XR,XR,-YR,YR);
inmesh=and(not(onedge),eq(sumfun2(wt,-1,1,-1,1),sumfun2(maxwt,-1,1,-1,1)));
lo=if(inmesh,1e10,if(wts-EPS,sum/wts,($lim[5]+$lim[6])/2));
_EOF_
# Generate height image
xform -m void $mesh | oconv - /tmp/edge$$.rad > /tmp/mesh$$.oct
cnt $res[2] $res[1] \
	| rcalc -e 'sp=$2/'"($res[1]-1)" -e 'tp=1-$1/'"($res[2]-1)" \
	-f /tmp/mesh$$.cal -e "xp=xpos(sp,tp);yp=ypos(sp,tp)" \
	-e '$1=xp;$2=yp;$3=z0;$4=0;$5=0;$6=1' \
	| rtrace -w -faf -oL -x $res[1] -y $res[2] /tmp/mesh$$.oct \
	| pvalue -r -b -df > /tmp/mesh$$.pic
rm /tmp/edge$$.rad /tmp/mesh$$.oct
# Output mesh surround
cat << _EOF_
# $0 $argv[*]
# Corner points:
#	$c00
#	$c01
#	$c11
#	$c10
_EOF_
pflip -v /tmp/mesh$$.pic \
	| pcomb -f /tmp/mesh$$.cal - \
	| pvalue -h -H -pG -d \
	| gensurf $mat $nam 'xpos(s,t)' 'ypos(s,t)' - \
		`ev $res[1]-1 $res[2]-1` \
		-e 'valid(xf,yf)=1e9-Z`SYS(xf,yf)' \
		-f /tmp/mesh$$.cal $smooth
rm /tmp/mesh$$.cal /tmp/mesh$$.pic
# All done -- exit
exit 0
# Usage error message
userr:
echo Usage: $0 "mat name mesh_file [-c corners][-b bord_width][-r step_size][-s]"
exit 1
