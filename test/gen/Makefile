# RCSid $Id: Makefile,v 1.2 2019/01/08 01:11:32 greg Exp $
#
# Unit tests for tools built in src/gen but not covered by test/renders
#

all:	test-gensky test-gendaymtx

clean:
	rm -f NYC11-15-14-21.rad test.smx

test-gensky:
	gensky 11 15 14:21EST +s -g .25 -t 3.5 \
-a 40.7128 -o 74.006 > NYC11-15-14-21.rad
	radcompare ref/NYC11-15-14-21.rad NYC11-15-14-21.rad

test-gendaymtx:	test.wea
	gendaymtx -r 90 -m 1 -g .3 .2 .1 -c .9 .9 1.2 test.wea > test.smx
	radcompare ref/test.smx test.smx
