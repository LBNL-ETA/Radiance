# RCSid $Id: Makefile,v 1.2 2019/01/26 00:55:50 greg Exp $
#
# Unit tests for tools built in src/util but not covered by test/renders
#

all:	test-vwright test-getinfo test-rcollate test-rmtxop test-dctimestep \
test-genskyvec

clean:
	rm -f test.mtx

test-vwright:	test.vf
	vwright -vf test.vf 3.5 > vwright.txt
	radcompare ref/vwright.txt vwright.txt
	rm -f vwright.txt

test-getinfo:	test.mtx
	getinfo -a Guppies "Fredo the Frog" < test.mtx | getinfo > getinfo.txt
	radcompare ref/getinfo.txt getinfo.txt
	rm -f getinfo.txt
	
test-rmtxop:	test.mtx
	rmtxop -ff -c .3 .9 .2 test.mtx -c .7 .2 .3 -t test.mtx > rmtxop.mtx
	radcompare ref/rmtxop.mtx rmtxop.mtx
	rm -f rmtxop.mtx

test-rcollate:	test.mtx
	radcompare ref/test.mtx test.mtx

test.mtx:	../renders/ref/rfmirror.mtx
	rcollate -t ../renders/ref/rfmirror.mtx > test.mtx

test-genskyvec:	test-dctimestep

test-dctimestep:	test.mtx
	gensky 3 21 10:15PST +s -g .3 -g 2.5 -a 36 -o 124 \
		| genskyvec -m 1 -c .92 1.03 1.2 \
		| dctimestep '!rmtxop -ff -t test.mtx' > dctimestep.mtx
	radcompare ref/dctimestep.mtx dctimestep.mtx
	rm -f dctimestep.mtx
