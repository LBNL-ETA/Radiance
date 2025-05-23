# RCSid $Id: Makefile,v 1.3 2025/05/23 19:42:48 greg Exp $
#
# Unit tests for tools built in src/cv
#

# Leaving off test-checkBSDF due to random sampling issues

all:	test-obj2rad test-robjutil test-mgf2rad test-rad2mgf

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
	radcompare -n -c# ref/sawtooth.mgf sawtooth.mgf 

test-mgf2rad:	testspecmat.mgf sawtooth.mgf
	mgf2rad -s testspecmat.mgf sawtooth.mgf > reconv.rad
	radcompare -c# ref/reconv.rad reconv.rad
	rm -f reconv.rad

sawtooth.mgf:	sawtooth.rad sawtooth.pts
	rad2mgf sawtooth.rad > sawtooth.mgf

rounded_box.obj:
	genbox plain box 5 8 3.5 -r .2 -n 6 -s -o > rounded_box.obj
