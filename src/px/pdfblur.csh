#!/bin/csh -f
# RCSid: $Id: pdfblur.csh,v 2.6 2005/01/18 01:37:43 greg Exp $
#
# Generate views for depth-of-field blurring on picture
#
if ($#argv != 3) then
	echo "Usage: $0 aperture nsamp viewfile"
	exit 1
endif
set a = "$1"
set n = "$2"
set vf = "$3"
cnt $n | rcalc -e `vwright i < $vf` \
-e "M:$n/5+1;a:$a/2;N:$n;" -e 'tmax:PI*a*(M+1)' \
-e 't=tmax/N*($1+rand($1))' \
-e 'theta=2*M*PI/(M-1)*(M-sqrt(M*M-(M-1)/(PI*a)*t))' \
-e 'r=a*(1-(M-1)/(2*M*M*PI)*theta)' \
-e 'rcost=r*cos(theta);rsint=r*sin(theta)' \
-e 'opx=ipx+rcost*ihx+rsint*ivx' \
-e 'opy=ipy+rcost*ihy+rsint*ivy' \
-e 'opz=ipz+rcost*ihz+rsint*ivz' \
-e 'os=is-rcost/(id*ihn);ol=il-rsint/(id*ivn)' \
-o 'VIEW= -vp ${opx} ${opy} ${opz} -vs ${os} -vl ${ol}'
