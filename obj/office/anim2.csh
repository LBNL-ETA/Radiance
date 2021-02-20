#!/bin/csh

set tempfile=/usr/tmp/pict.$$
set destdir=/usr2/greg/ras/anim/model2
set nframes=$1
mkdir $destdir
set i=0
while ($i < $nframes)
	set view=(`ev $i/$nframes\*6.28 | rcalc -f anim2.cal`)
	rpict -vp $view[1-3] -vd $view[4-6] -vu 0 1 0 -vh 38 -vv 30 \
		-x 512 -y 400 -t 3600 -e anim2.err model.oct \
	| pfilt -1 -x 512 -y 400 -e -2 > $tempfile
	ra_pr $tempfile $destdir/$i
	rm $tempfile
	@ i++
end
