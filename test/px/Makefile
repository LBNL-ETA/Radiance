# RCSid $Id: Makefile,v 1.7 2020/07/26 20:15:33 greg Exp $
#
# Unit tests for tools built in src/px
#

HDR_CMP = radcompare -rms 0.07 -rel 1e-6 -max 0.2

IMG_CMP = radcompare -h -rms 0.03 -rel 0.007 -max 0.2

all:	test-pfilt test-pvalue test-ra_tiff test-pcompos \
test-protate test-pflip test-ra_bmp test-ra_xyze test-pcond test-ra_t16 \
test-ra_rgbe test-ra_ppm test-pextrem test-pcomb test-pinterp test-psign \
test-falsecolor test-ra_ps

clean:
	rm -f test.hdr test24.hdr test24.tif test32L.tif pcomb.hdr \
norm24.tif norm8.tif pcompos.hdr pcondh.hdr test_extrema.txt psign.hdr \
pinterp.hdr false.hdr ra_psC.ps

test-pcond:	test.hdr
	pcond -h test.hdr pcondh.hdr
	$(IMG_CMP) ref/pcondh.hdr pcondh.hdr
	rm -f pcondh.hdr

test-pcomb:	test.hdr
	pcomb -e 'm=1/(gi(1)+.1);ro=gi(1)*m;go=bi(1)*m;bo=ri(1)*m' \
test.hdr > pcomb.hdr
	$(HDR_CMP) ref/pcomb.hdr pcomb.hdr
	rm -f pcomb.hdr

test-ra_xyze:	test-pfilt

test-pfilt:	test.hdr
	$(HDR_CMP) ref/test.hdr test.hdr

test-pvalue:	test24.hdr
	$(IMG_CMP) ref/test24.hdr test24.hdr

test-pextrem:	test.hdr
	pextrem -o test.hdr > test_extrema.txt
	radcompare ref/test_extrema.txt test_extrema.txt
	rm -f test_extrema.txt

test-protate test-pflip:	test-pcompos

test-pcompos:	test.hdr test24.hdr
	protate -r test.hdr | pcompos -a 2 -s 10 -b 1 1 1 \
		'!pflip -h test24.hdr' - \
		'!pflip -v test24.hdr' '!protate test.hdr' > pcompos.hdr
	$(HDR_CMP) ref/pcompos.hdr pcompos.hdr
	rm -f pcompos.hdr

test-pinterp:	test.hdr
	pinterp -vtv -vh 60 -vv 40 -x 200 -y 200 test.hdr 1.0 > pinterp.hdr
	$(HDR_CMP) ref/pinterp.hdr pinterp.hdr
	rm -f pinterp.hdr

test-psign:
	psign -cb .9 .4 .01 -cf .1 .2 1 -dd -h 200 -a 1.5 -s -.1 Testing > psign.hdr
	$(HDR_CMP) ref/psign.hdr psign.hdr
	rm -f psign.hdr

test-falsecolor:	test.hdr
	perl ../../src/px/falsecolor.pl -s 2e5 -log 5 -i test.hdr > false.hdr
	$(HDR_CMP) ref/false.hdr false.hdr
	rm -f false.hdr

test-ra_rgbe:	test.hdr
	ra_rgbe -e -5 test.hdr | ra_rgbe -r -e +5 | $(HDR_CMP) -w test.hdr -

test-ra_tiff:	test.hdr test24.hdr
	ra_tiff test24.hdr test24.tif
	ra_tiff -r test24.tif | $(IMG_CMP) test24.hdr -
	ra_tiff -L test.hdr test32L.tif
	ra_tiff -r test32L.tif | $(HDR_CMP) -h test.hdr -
	rm -f test24.tif test32L.tif

test-normtiff:	test.hdr
	normtiff -z test.hdr norm24.tif
	radcompare ref/norm24.tif norm24.tif
	normtiff -z -b test.hdr norm8.tif
	radcompare ref/norm8.tif norm8.tif
	rm -f norm24.tif norm8.tif

test-ra_bmp:	test24.hdr
	ra_bmp test24.hdr | ra_bmp -r | $(IMG_CMP) test24.hdr -

test-ra_t16:	test24.hdr
	ra_t16 -3 test24.hdr | ra_t16 -r | $(IMG_CMP) test24.hdr -

test-ra_ppm:	test24.hdr
	ra_ppm test24.hdr | ra_ppm -r | $(IMG_CMP) test24.hdr -

test-ra_ps:	test.hdr
	ra_ps -b -C test.hdr > ra_psC.ps
	radcompare ref/ra_psC.ps ra_psC.ps
	rm -f ra_psC.ps

test24.hdr:	test.hdr
	pvalue -g 2.2 -db test.hdr | pvalue -r -g 2.2 -db > test24.hdr

test.hdr:	../renders/sunset.hdr
	pfilt -x 544 -y 544 -r .6 ../renders/sunset.hdr \
| ra_xyze -r > test.hdr
