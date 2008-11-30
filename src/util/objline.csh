#!/bin/csh -f
# RCSid: $Id: objline.csh,v 2.6 2008/11/30 19:56:02 greg Exp $
# Create four standard views of scene and present as line drawings
#
set oblqxf="-rz 45 -ry 45"
onintr quit
set d=`mktemp -d /tmp/ol.XXXXXX`
if ($#argv) then
	set origf=""
	set oblqf=""
	foreach f ($argv)
		set fn="$f:t"
		set fn="$fn:r"
		xform $f:q > $d/$fn:q.orig
		rad2mgf $d/$fn:q.orig > $d/$fn:q.orig.mgf
		set origf=($origf:q $fn:q.orig.mgf)
		echo i $fn:q.orig.mgf $oblqxf > $d/$fn:q.oblq.mgf
		set oblqf=($oblqf:q $fn:q.oblq.mgf)
	end
else
	set origf=stdin.orig.mgf
	set oblqf=stdin.oblq.mgf
	xform > $d/stdin.orig
	rad2mgf $d/stdin.orig > $d/stdin.orig.mgf
	echo i stdin.orig.mgf $oblqxf > $d/stdin.oblq.mgf
endif
cd $d
set rce='xm=($1+$2)/2;ym=($3+$4)/2;zm=($5+$6)/2;\
max(a,b):if(a-b,a,b);r=max(max($2-$1,$4-$3),$6-$5)*.52;\
$1=xm-r;$2=xm+r;$3=ym-r;$4=ym+r;$5=zm-r;$6=zm+r'
set origdim=`getbbox -h *.orig | rcalc -e $rce:q`
set oblqdim=`xform $oblqxf *.orig | getbbox -h | rcalc -e $rce:q`
mgf2meta -t .005 x $origdim $origf > x.mta
mgf2meta -t .005 y $origdim $origf > y.mta
mgf2meta -t .005 z $origdim $origf > z.mta
mgf2meta -t .005 x $oblqdim $oblqf > o.mta
plot4 {x,y,z,o}.mta
quit:
cd
exec rm -rf $d
