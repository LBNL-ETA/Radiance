# RCSid $Id: Makefile,v 1.3 2019/01/08 23:53:04 greg Exp $
#
# Unit tests for tools built in src/gen but not covered by test/renders
#

all:	test-gensky test-gendaymtx test-genbox test-genrev test-genworm \
test-genprism test-gensurf

clean:
	rm -f gen*.rad

test-gensky:
	gensky 11 15 14:21EST +s -g .25 -t 3.5 \
-a 40.7128 -o 74.006 > NYC11-15-14-21.rad
	radcompare ref/NYC11-15-14-21.rad NYC11-15-14-21.rad
	rm -f NYC11-15-14-21.rad

test-gendaymtx:	test.wea
	gendaymtx -r 90 -m 1 -g .3 .2 .1 -c .9 .9 1.2 test.wea > test.smx
	radcompare ref/test.smx test.smx
	rm -f test.smx

test-genbox:
	genbox tmat tbox 1 2 3 -i -b .1 > genbox.rad
	radcompare ref/genbox.rad genbox.rad
	rm -f genbox.rad

test-genrev:
	genrev tmat trev 'sin(2*PI*t)' '2+cos(2*PI*t)' 16 -s > genrev.rad
	radcompare ref/genrev.rad genrev.rad
	rm -f genrev.rad

test-genworm:
	genworm tmat tworm '0' '5*sin(t)' '5*cos(t)' '.4-(.5-t)*(.5-t)' 8 > genworm.rad
	radcompare ref/genworm.rad genworm.rad
	rm -f genworm.rad

test-genprism:
	genprism tmat tprism 8 0 5 4 5 4 4 24.5 4 24.5 3 30 1.5 30 22 0 22 \
-l 0 0 -1.5 -r .2 > genprism.rad
	radcompare ref/genprism.rad genprism.rad
	rm -f genprism.rad

test-gensurf:
	gensurf tmat tsurf '15.5+x(theta(s),phi(t))' '10.5+y(theta(s),phi(t))' \
'30.75+z(theta(s),phi(t))' 4 15 -f basin.cal -s > gensurf.rad
	radcompare ref/gensurf.rad gensurf.rad
	rm -f gensurf.rad

