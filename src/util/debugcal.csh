#!/bin/csh -f
# RCSid: $Id: debugcal.csh,v 2.8 2021/08/26 18:00:39 greg Exp $
#
# Script to debug cal files for Radiance
#
# Takes octree and rcalc options as arguments.
# One of the rcalc options should be -f calfile.
# Note that the real arguments arg(n) or A1, A2, etc must also be given.
# Input is a ray origin and direction, such as that produced by ximage.
#
if ( $#argv < 2 ) then
	echo "Usage: $0 octree [rcalc options]"
	exit 1
endif

rtrace -h- -fad -x 1 -odNplLcn $1 | rcalc -id16 -u \
		-e 'Dx=$1;Dy=$2;Dz=$3;Nx=$4;Ny=$5;Nz=$6;Px=$7;Py=$8;Pz=$9' \
		-e 'T=$10;Ts=$11;Lu=$12;Lv=$13;NxP=$14;NyP=$15;NzP=$16'\
		-e 'S:1;Tx:0;Ty:0;Tz:0:Ix:1;Iy:0;Iz:0;Jx:0;Jy:1;Jz:0;Kx:0;Ky:0;Kz:1' \
		-e 'Rdot=-Dx*Nx-Dy*Ny-Dz*Nz' -e 'RdotP=-Dx*NxP-Dy*NyP-Dz*NzP' \
		-e 'CrP=A1;CgP=A2;CbP=A3' -e 'DxA:0;DyA:0;DzA:0' \
		-f rayinit.cal $argv[2-]:q
