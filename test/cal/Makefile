# RCSid $Id: Makefile,v 1.2 2020/06/29 22:46:22 greg Exp $
#
# Tests for tools built in ray/src/cal directory.
#

all:	test-cnt test-rcalc test-total test-histo test-rlam \
test-tabfunc test-rsplit test-ev

clean:
	rm -f cnt.txt rcalc.txt total.txt histo.txt rlam.txt \
tabfunc.cal tabfunc.txt cfirst.txt cthird.txt evout.txt

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

test-tabfunc:	tabfunc.txt
	radcompare ref/tabfunc.txt tabfunc.txt

test-rsplit:	cfirst.txt cthird.txt
	radcompare ref/cfirst.txt cfirst.txt
	radcompare ref/cthird.txt cthird.txt

test-ev:	evout.txt
	radcompare ref/evout.txt evout.txt

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

tabfunc.txt:	tabfunc.cal
	cnt  50 | rcalc -f tabfunc.cal -e 'xv=$$1+5' \
-e '$$1=cfirst(xv);$$2=csecond(xv);$$3=cthird(xv)' > tabfunc.txt

tabfunc.cal:	histo.txt
	tabfunc -i cfirst csecond cthird < histo.txt > tabfunc.cal

cfirst.txt cthird.txt:	tabfunc.txt
	rsplit -on 11 -f '-t ' - . cthird.txt < tabfunc.txt > cfirst.txt

evout.txt:
	ev 'sin(0.9353)' 'sqrt(2)' 'log(atan(.5)+3.917)' > evout.txt
