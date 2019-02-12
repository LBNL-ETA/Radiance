# RCSid $Id: Makefile,v 1.1 2019/02/12 22:54:19 greg Exp $
#
# Tests for tools built in ray/src/cal directory.
#

all:	test-cnt test-rcalc test-total test-histo test-rlam

clean:
	rm -f cnt.txt rcalc.txt total.txt histo.txt rlam.txt

test-cnt:	cnt.txt
	radcompare ref/cnt.txt cnt.txt

test-rcalc:	rcalc.txt
	radcompare ref/rcalc.txt rcalc.txt

test-total:	total.txt
	radcompare ref/total.txt total.txt

test-histo:	histo.txt
	radcompare ref/histo.txt histo.txt

test-rlam:	rlam.txt
	radcompare ref/rlam.txt rlam.txt

cnt.txt:
	cnt 5 3 2 > cnt.txt

rcalc.txt:	cnt.txt
	rcalc -o 'Test $${v1} $$(s1) $${v2}' \
-e 'v1=$$1*$$2;v2=($$2-$$1)*exp($$3)' -s s1=HEY cnt.txt > rcalc.txt

total.txt:	cnt.txt
	total cnt.txt > total.txt
	total -l cnt.txt >> total.txt
	total -u cnt.txt >> total.txt
	total -m cnt.txt >> total.txt
	total -s2.5 cnt.txt >> total.txt
	total -3 -r cnt.txt >> total.txt

histo.txt:	total.txt
	histo 0 60 5 < total.txt > histo.txt

rlam.txt:	total.txt cnt.txt histo.txt
	rlam -in 5 total.txt cnt.txt histo.txt > rlam.txt
