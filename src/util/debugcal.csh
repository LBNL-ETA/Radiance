#!/bin/csh -f
# RCSid: $Id: debugcal.csh,v 2.6 2003/02/22 02:07:30 greg Exp $
#
# Script to debug cal files for Radiance
#
# Takes octree and rcalc options as arguments.
# One of the rcalc options should be -f calfile.
# Note that the real arguments A1, A2, etc must also be given.
# Input is a ray origin and direction, such as that produced by ximage.
#
if ( $#argv < 2 ) then
	echo "Usage: $0 octree [rcalc options]"
	exit 1
endif
if ( ! $?RAYPATH ) then
	set RAYPATH=.:/usr/local/lib/ray
endif
set initfile=
foreach d (`echo $RAYPATH | sed 's/:/ /g'`)
	if ( -r $d/rayinit.cal ) then
		set initfile="-f $d/rayinit.cal"
		break
	endif
end

rtrace -h- -x 1 -odNplL $1 | rcalc -u -e 'Dx=$1;Dy=$2;Dz=$3' \
		-e 'Nx=$4;Ny=$5;Nz=$6;Px=$7;Py=$8;Pz=$9' \
		-e 'T=$10;Ts=$11' -e 'S:1;Tx:0;Ty:0;Tz:0' \
		-e 'Ix:1;Iy:0;Iz:0;Jx:0;Jy:1;Jz:0;Kx:0;Ky:0;Kz:1' \
		-e 'Rdot=-Dx*Nx-Dy*Ny-Dz*Nz' -e 'RdotP=Rdot' \
		-e 'NxP=Nx;NyP=Ny;NzP=Nz' -e 'CrP=A1;CgP=A2;CbP=A3' \
		-e 'DxA:0;DyA:0;DzA:0' \
		$initfile $argv[2-]:q
