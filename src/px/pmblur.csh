#!/bin/csh -f
# RCSid: $Id: pmblur.csh,v 2.2 2003/02/22 02:07:27 greg Exp $
#
# Generate views for motion blurring on picture
#
if ($#argv != 4) then
	echo "Usage: $0 speed nsamp v0file v1file" >/dev/tty
	exit 1
endif
set s = "$1"
set n = "$2"
set vc = "$3"
set vn = "$4"
cnt $n | rcalc -e `vwright C < $vc` -e `vwright N < $vn` \
-e "t=$s/$n"'*($1+rand($1))' \
-e 'opx=(1-t)*Cpx+t*Npx;opy=(1-t)*Cpy+t*Npy;opz=(1-t)*Cpz+t*Npz' \
-e 'odx=(1-t)*Cdx+t*Ndx;ody=(1-t)*Cdy+t*Ndy;odz=(1-t)*Cdz+t*Ndz' \
-e 'oux=(1-t)*Cux+t*Nux;ouy=(1-t)*Cuy+t*Nuy;ouz=(1-t)*Cuz+t*Nuz' \
-e 'oh=(1-t)*Ch+t*Nh;ov=(1-t)*Cv+t*Nv' \
-o 'VIEW= -vp ${opx} ${opy} ${opz} -vd ${odx} ${ody} ${odz} -vu ${oux} ${ouy} ${ouz} -vh ${oh} -vv ${ov}'
