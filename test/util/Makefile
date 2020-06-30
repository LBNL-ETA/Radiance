# RCSid $Id: Makefile,v 1.6 2020/06/30 18:00:19 greg Exp $
#
# Unit tests for tools built in src/util but not covered by test/renders
#

all:	test-vwright test-getinfo test-rcollate test-rmtxop test-dctimestep \
test-genskyvec test-vwrays test-glarendx test-rcode_norm test-rcode_depth \
test-rcode_ident test-wrapBSDF test-rttree_reduce

clean:
	rm -f test.mtx ttree.txt vwrays.txt

test-vwright:	test.vf
	vwright -vf test.vf 3.5 > vwright.txt
	radcompare ref/vwright.txt vwright.txt
	rm -f vwright.txt

test-vwrays:	vwrays.txt
	radcompare ref/vwrays.txt vwrays.txt

test-getinfo:	test.mtx
	getinfo -a Guppies "Fredo the Frog" < test.mtx | getinfo > getinfo.txt
	radcompare ref/getinfo.txt getinfo.txt
	rm -f getinfo.txt
	
test-rmtxop:	test.mtx
	rmtxop -ff -c .3 .9 .2 test.mtx -c .7 .2 .3 -t test.mtx | radcompare ref/rmtxop.mtx - 

test-rcollate:	test.mtx
	radcompare ref/test.mtx test.mtx

test-rcode_ident:	idents.txt
	rcode_ident -8 -x 20 -y 10 idents.txt | radcompare ref/idents.idx -
	rcode_ident -r -H ref/idents.idx | radcompare idents.txt -

test-rcode_norm:	vwrays.txt
	rcode_norm -hi -x 40 -y 15 vwrays.txt | radcompare ref/norms.nrm -

test-rcode_depth:	depths.txt
	rcode_depth -x 20 -y 10 -d 150/inch depths.txt | radcompare -h ref/depths.dpt -
	rcode_depth -r -Ho ref/depths.dpt | radcompare -w depths.txt -

test-glarendx:	glaresrc.txt
	glarendx -t guth_vcp glaresrc.txt | radcompare ref/guth_vcp.txt -
	glarendx -t cie_cgi glaresrc.txt | radcompare ref/cie_cgi.txt -
	glarendx -t ugr glaresrc.txt | radcompare ref/ugr.txt -
	glarendx -t brs_gi glaresrc.txt | radcompare ref/brs_gi.txt -
	glarendx -t dgi glaresrc.txt | radcompare ref/dgi.txt -
	glarendx -t guth_dgr glaresrc.txt | radcompare ref/guth_dgr.txt -
	glarendx -t vert_dir glaresrc.txt | radcompare ref/vert_dir.txt -
	glarendx -t vert_ill glaresrc.txt | radcompare ref/vert_ill.txt -
	glarendx -t vert_ind glaresrc.txt | radcompare ref/vert_ind.txt -

test-rttree_reduce: ttree.txt
	radcompare ref/ttree.txt ttree.txt

test-wrapBSDF:	ttree.txt
	wrapBSDF -a t3 -u centimeter -f 'm=Test Guys;n=Bogus;t=2.54;h=20;w=45;' \
-rb ttree.txt > testBSDF.xml
	radcompare ref/testBSDF.xml testBSDF.xml
	rm -f testBSDF.xml

test.mtx:	../renders/ref/rfmirror.mtx
	rcollate -t ../renders/ref/rfmirror.mtx > test.mtx

ttree.txt:
	cnt 256 | rcalc -e '$$1=floor(sin(PI/256*recno)*40)' \
| rttree_reduce -h -g 3 -r 3 -t 95 > ttree.txt

vwrays.txt:
	vwrays -vts -vp 100 50 -75 -vo 10 -va 1000 -vh 160 -vv 130 \
-x 20 -y 20 -vu 0 -1 0 -vd .358 -.8851 1.152 > vwrays.txt

test-genskyvec:	test-dctimestep

test-dctimestep:	test.mtx
	gensky 3 21 10:15PST +s -g .3 -g 2.5 -a 36 -o 124 \
		| genskyvec -m 1 -c .92 1.03 1.2 \
		| dctimestep '!rmtxop -ff -t test.mtx' > dctimestep.mtx
	radcompare ref/dctimestep.mtx dctimestep.mtx
	rm -f dctimestep.mtx
