#!/bin/csh -f
# RCSid: $Id: objline.csh,v 2.5 2008/08/25 04:50:32 greg Exp $
# Create four standard views of scene and present as line drawings
#
set oblqxf="-rz 45 -ry 45"
onintr quit
set d=`mktemp -d /tmp/ol.XXXXXX`
if ($#argv) then
	set origf=""
	set oblqf=""
	foreach f ($argv)
		xform -e $f > $d/$f.orig
		rad2mgf $d/$f.orig > $d/$f:r.orig.mgf
		set origf=($origf $f:r.orig.mgf)
		echo i $f:r.orig.mgf $oblqxf > $d/$f:r.oblq.mgf
		set oblqf=($oblqf $f:r.oblq.mgf)
	end
else
	set origf=stdin.orig.mgf
	set oblqf=stdin.oblq.mgf
	xform -e > $d/stdin.orig
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
