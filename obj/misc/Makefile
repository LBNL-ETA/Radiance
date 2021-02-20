#
#  Makefile for examples
#

DEV=x11

daf:	daf.oct
	rview -o $(DEV) -vp 50 60 40 -vd 0 -0.707107 -0.707107 -vu 0 1 0 \
-vh 20 -vv 20 -av .2 .2 .2 daf.oct

daffodil.hdr:	daf.oct
	rpict -vp 50 60 40 -vd 0 -0.707107 -0.707107 -vu 0 1 0 \
-vh 20 -vv 20 -x 470 -y 470 -av .2 .2 .2 daf.oct > daffodil.hdr

daf.oct:	daffodil.xf room
	oconv daffodil.xf room > daf.oct

daffodil.xf:	daffodil
	xform -e -rz -40 -t 50 40 20 daffodil > daffodil.xf

exm:	examp.oct
	rview -o $(DEV) -vp 3 2 5 -vd -0.268328 -0.357771 -0.894427 \
-vu 0 1 0 -vh 45 -vv 45 -av .25 .25 .25 examp.oct

examp.hdr:	examp.oct
	rpict -vp 3 2 5 -vd -0.268328 -0.357771 -0.894427 \
-vu 0 1 0 -vh 45 -vv 45 -x 470 -y 470 -av .5 .5 .5 \
-ab 0 -ad 512 -as 128 examp.oct > examp.hdr

examp.oct:	examp
	oconv examp > examp.oct
