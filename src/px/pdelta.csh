#!/bin/csh -f
# RCSid: $Id: pdelta.csh,v 3.3 2008/11/10 19:08:19 greg Exp $
#
# Compute 1976 CIE Lab deltaE* between two Radiance pictures
#
if ($#argv != 2) then
	echo "Usage: $0 picture1 picture2 > output.hdr"
	exit 1
endif
set cielab='sq(x):x*x;Ls(Yi):if(Yi/Yw-.01,116*(Yi/Yw)^(1/3)-16,903.3*Yi/Yw);as(Xi,Yi,Zi):500*((Xi/Xw)^(1/3)-(Yi/Yw)^(1/3));bs(Xi,Yi,Zi):200*((Yi/Yw)^(1/3)-(Zi/Zw)^(1/3));dE(X1,Y1,Z1,X2,Y2,Z2):sqrt(sq(Ls(Y1)-Ls(Y2))+sq(as(X1,Y1,Z1)-as(X2,Y2,Z2))+sq(bs(X1,Y1,Z1)-bs(X2,Y2,Z2)));'
# The following is for Radiance RGB -> XYZ
set rgb2xyz='X(R,G,B):92.03*R+57.98*G+28.99*B;Y(R,G,B):47.45*R+119.9*G+11.6*B;Z(R,G,B):4.31*R+21.99*G+152.7*B'
# The following is for sRGB -> XYZ
# set rgb2xyz='X(R,G,B):73.82*R+64.01*G+32.31*B;Y(R,G,B):38.06*R+128*G+12.92*B;Z(R,G,B):3.46*R+21.34*G+170.15*B'
set f1="$1"
set inp1='x1=ri(1);y1=gi(1);z1=bi(1)'
set f2="$2"
set inp2='x2=ri(2);y2=gi(2);z2=bi(2)'
if ( "`getinfo < $f1:q | grep '^FORMAT=32-bit_rle_xyze'`" == "" ) then
	set inp1='x1=X(ri(1),gi(1),bi(1));y1=Y(ri(1),gi(1),bi(1));z1=Z(ri(1),gi(1),bi(1))'
endif
if ( "`getinfo < $f2:q | grep '^FORMAT=32-bit_rle_xyze'`" == "" ) then
	set inp2='x2=X(ri(2),gi(2),bi(2));y2=Y(ri(2),gi(2),bi(2));z2=Z(ri(2),gi(2),bi(2))'
endif
set tempf=`mktemp /tmp/tf.XXXXXX`
pfilt -1 -x 128 -y 128 -p 1 $f1:q | pvalue -o -h -H -d > $tempf
set wht=(`total -u $tempf`)
set avg=`rcalc -e '$1=$2' $tempf | total -m`
rm $tempf
pcomb -e $cielab:q -e $rgb2xyz:q \
	-e "Yw:179*3*$avg; Xw:$wht[1]*Yw/$wht[2]; Zw:$wht[3]*Yw/$wht[2]" \
	-e $inp1:q -e $inp2:q -e 'lo=dE(x1,y1,z1,x2,y2,z2)' \
	-o $f1:q -o $f2:q
