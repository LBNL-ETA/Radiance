#!/bin/csh -f
# RCSid $Id: pmdblur.csh,v 3.1 2005/01/18 03:59:41 greg Exp $
#
# Generate views for motion and depth blurring on picture
#
if ($#argv != 5) then
	echo "Usage: $0 speed aperture nsamp v0file v1file"
	exit 1
endif
set s = "$1"
set a = "$2"
set n = "$3"
set vc = "$4"
set vn = "$5"
if (`ev "if($s-.01,0,1)"`) then
	pdfblur $a $n $vc
	exit
endif
if (`ev "if($a,0,1)"`) then
	pmblur $s $n $vc $vn
	exit
endif
cnt $n | rcalc -e `vwright C < $vc` -e `vwright N < $vn` \
-e "t=$s/$n"'*($1+rand($1))' \
-e "r=$a/2"'*sqrt(rand(182+7*$1));theta=2*PI*rand(-10-$1)' \
-e 'rcost=r*cos(theta);rsint=r*sin(theta)' \
-e 'opx= (1-t)*(Cpx+rcost*Chx+rsint*Cvx) + t*(Npx+rcost*Nhx+rsint*Nvx)' \
-e 'opy= (1-t)*(Cpy+rcost*Chy+rsint*Cvy) + t*(Npy+rcost*Nhy+rsint*Nvy)' \
-e 'opz= (1-t)*(Cpz+rcost*Chz+rsint*Cvz) + t*(Npz+rcost*Nhz+rsint*Nvz)' \
-e 'odx= (1-t)*Cdx*Cd + t*Ndx*Nd' \
-e 'ody= (1-t)*Cdy*Cd + t*Ndy*Nd' \
-e 'odz= (1-t)*Cdz*Cd + t*Ndz*Nd' \
-e 'oux=(1-t)*Cux+t*Nux;ouy=(1-t)*Cuy+t*Nuy;ouz=(1-t)*Cuz+t*Nuz' \
-e 'oh=(1-t)*Ch+t*Nh;ov=(1-t)*Cv+t*Nv' \
-e 'os= (1-t)*(Cs-rcost/(Cd*Chn)) + t*(Ns-rcost/(Nd*Nhn))' \
-e 'ol= (1-t)*(Cl-rsint/(Cd*Cvn)) + t*(Nl-rsint/(Nd*Nvn))' \
-o 'VIEW= -vp ${opx} ${opy} ${opz} -vd ${odx} ${ody} ${odz} -vu ${oux} ${ouy} ${ouz} -vh ${oh} -vv ${ov} -vs ${os} -vl ${ol}'
