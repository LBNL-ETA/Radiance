# RCSid $Id: Makefile,v 1.1 2020/06/30 23:52:50 greg Exp $
#
# Unit tests for tools built in src/px
#

HDR_CMP = radcompare -rms 0.07 -rel 1e-6 -max 0.2

IMG_CMP = radcompare -h -rms 0.03 -rel 0.03 -max 0.07

all:	test-pfilt test-pvalue test-ra_tiff test-normtiff test-pcompos \
test-protate test-pflip

clean:
	rm -f test.hdr test24.hdr test24.tif test32L.tif \
norm24.tif norm8.tif pcompos.hdr

test-pfilt:	test.hdr
	$(HDR_CMP) ref/test.hdr test.hdr

test-pvalue:	test24.hdr
	$(IMG_CMP) ref/test24.hdr test24.hdr

test-protate test-pflip test-pcompos:	test.hdr test24.hdr
	protate -r test.hdr | pcompos -a 2 -s 10 -b 1 1 1 \
		'!pflip -h test24.hdr' - \
		'!pflip -v test24.hdr' '!protate test.hdr' > pcompos.hdr
	$(HDR_CMP) ref/pcompos.hdr pcompos.hdr
	rm -f pcompos.hdr

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

test24.hdr:	test.hdr
	pvalue -db test.hdr | pvalue -r -db > test24.hdr

test.hdr:	../renders/sunset.hdr
	pfilt -x /2.5 -y /2.5 -r .6 ../renders/sunset.hdr > test.hdr
