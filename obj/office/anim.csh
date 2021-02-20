#!/bin/csh

set tempfile=/usr/tmp/pict.$$
set destfile=~/ray/spool/anim/model
set sunrise=6
set sunset=18
set nframes=12
set i=0
while ($i < $nframes)
	oconv -i model.dark.oct "\!gensky 6 17 `ev $sunrise+\($sunset-$sunrise\)\*$i/$nframes` | xform -rz -90 -rx -90" skywindow \
	| rpict -vp 8 36 -27 -vd -.56 -.23 .79 -vu 0 1 0 -vh 39.9 -vv 27.5 \
		-x 300 -y 203 -ps 1 -ab 1 -ad 32 -ar 64 -ds .5 \
		-t 3600 -e model.day.err \
	| pfilt -1 -x 300 -y 203 -e 1.5 > $tempfile
	ra_pr $tempfile $destfile.$i
	rm $tempfile
	@ i++
end
