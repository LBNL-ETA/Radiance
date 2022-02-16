# RCSid $Id: Makefile,v 1.1 2022/02/16 18:54:13 greg Exp $
#
# Unit tests for tools built in src/cv
#

all:	test-checkBSDF test-obj2rad test-robjutil \
test-mgf2rad test-rad2mgf

clean:
	rm -f checkBSDFout.txt rounded_box_o2r.rad rounded_box_ru.rad \
sawtooth.mgf rounded_box.obj

test-checkBSDF:	../renders/blinds20c.xml ../renders/sawtooth.xml
	checkBSDF ../renders/blinds20c.xml ../renders/sawtooth.xml \
> checkBSDFout.txt
	radcompare ref/checkBSDFout.txt checkBSDFout.txt
	rm -f checkBSDFout.txt

test-obj2rad:	rounded_box.obj
	obj2rad rounded_box.obj > rounded_box_o2r.rad
	radcompare -c# ref/rounded_box_o2r.rad rounded_box_o2r.rad
	rm -f rounded_box_o2r.rad

test-robjutil:	rounded_box.obj
	robjutil -n +r +T -x "-ry 10 -t -3 -2 0.5" rounded_box.obj \
> rounded_box_ru.rad
	radcompare -c# ref/rounded_box_ru.rad rounded_box_ru.rad
	rm -f rounded_box_ru.rad

test-rad2mgf:	sawtooth.mgf
	radcompare -c# ref/sawtooth.mgf sawtooth.mgf 

test-mgf2rad:	sawtooth.mgf
	mgf2rad sawtooth.mgf > sawtooth_rec.rad
	radcompare -c# ref/sawtooth_rec.rad sawtooth_rec.rad
	rm -f sawtooth_rec.rad

sawtooth.mgf:	sawtooth.rad sawtooth.pts
	rad2mgf sawtooth.rad > sawtooth.mgf

rounded_box.obj:
	genbox plain box 5 8 3.5 -r .2 -n 6 -s -o > rounded_box.obj
