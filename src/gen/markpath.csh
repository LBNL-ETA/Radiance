#!/bin/csh -f
# RCSid: $Id: markpath.csh,v 2.1 2003/02/22 02:07:23 greg Exp $
#
# Put right trianglar markers down a path at the given intervals.
# Use with replmarks to place regular-sized objects along a path.
# Shorter (y-side) of triangles is always horizontal (perp. to z).
#
# Input is an ordered list of 3-D points defining the path.
# We interpolate the path and align our markers with it.
# Triangles are sized and aligned so tip of one meets butt of next
#
if ($#argv < 2) then
	echo "Usage: $0 3d.pts spacing [markmat]"
	exit 1
endif
set pts=$1
set step=$2
set mat=mark
if ($#argv > 2) set mat=$3
set npts=`wc -l < $pts`
(head -1 $pts ; cat $pts) | lam - $pts \
	| rcalc -e '$1=d($1,$2,$3,$4,$5,$6)' -e "cond=$npts+.5-recno" \
	-e 's(x):x*x;d(x0,y0,z0,x1,y1,z1):sqrt(s(x1-x0)+s(y1-y0)+s(z1-z0))' \
	| total -1 -r > /tmp/run$$.dat
lam /tmp/run$$.dat $pts | tabfunc -i xp yp zp > /tmp/path$$.cal
set tmax=`tail -1 /tmp/run$$.dat`
set nsteps=`ev "floor($tmax/$step)"`
echo $mat > /tmp/tri$$.fmt
cat >> /tmp/tri$$.fmt << '_EOF_'
polygon marker.${recno}
0
0
9
	${   x0    } ${   y0    } ${   z0   }
	${   x1    } ${   y1    } ${   z1   }
	${   x2    } ${   y2    } ${   z2   }

'_EOF_'
cnt $nsteps | rcalc -o /tmp/tri$$.fmt -f /tmp/path$$.cal -e st=$step \
	-e 't=$1*st;x0=xp(t);y0=yp(t);z0=zp(t)' \
	-e 'x1=xp(t+st);y1=yp(t+st);z1=zp(t+st)' \
	-e 'x2=x0+.5*(y0-y1);y2=y0+.5*(x1-x0);z2=z0'
rm /tmp/run$$.dat /tmp/path$$.cal /tmp/tri$$.fmt
