#!/bin/csh -f
# RCSid: $Id: genpine.csh,v 2.3 2005/02/16 05:40:11 greg Exp $
#
# Generate a tree
# Pine version 2
#
# First send header and parse arguments
#
onintr done
echo \# $0 $*
set nleaves=150
set nlevels=4
set aspect=1.2
unset needles
while ($#argv > 0)
	switch ($argv[1])
	case -r:
		shift argv
		set nlevels=$argv[1]
		breaksw
	case -n:
		shift argv
		set nleaves=$argv[1]
		breaksw
	case -o:
		shift argv
		set needles=$argv[1]
		breaksw
	case -a:
		shift argv
		set aspect=$argv[1]
		breaksw
	default:
		echo bad option $argv[1]
		exit 1
	endsw
	shift argv
end
#
#  Send materials
#
cat << _EOF_

void plastic bark_mat
0
0
5 .6 .5 .45 0 0

void plastic leaf_mat
0
0
5 .11 .36 .025 0 0
_EOF_
#
# Next start seedling
#
set tree=/tmp/t$$
set oldtree=/tmp/ot$$
set thisrad=.035
cat << _EOF_ > $tree

void colorpict bark_pat
9 red green blue pinebark.pic cyl.cal cyl_match_u cyl_match_v -s $thisrad
0
2 1.5225225 1

bark_pat alias my_bark_mat bark_mat

my_bark_mat cone top
0
0
8
	0	0	0
	0	0	1
	$thisrad	.02

my_bark_mat sphere tip
0
0
4	0	0	1	.02
_EOF_
if ( ! $?needles ) set needles=n.$nleaves.oct
if ( ! -f $needles ) then
oconv -f "\!cnt $nleaves | rcalc -e nl=$nleaves -o needle.fmt" > $needles
endif
echo leaf_mat instance needles 1 $needles 0 0 >> $tree
#
# Now grow tree:
#
#	1) Save oldtree
#	2) Move tree up and extend trunk
#	3) Duplicate oldtree at branch positions
#	4) Repeat
#
@ i=0
while ($i < $nlevels)
	mv -f $tree $oldtree
	set lastrad=$thisrad
	set move=`ev "(2*$aspect)^($i+1)"`
	set thisrad=`ev "$lastrad+$move*.015"`
	xform -ry `ev "25/($i+1)"` -t 0 0 $move $oldtree > $tree
	echo void colorpict bark_pat 9 red green blue pinebark.pic \
			cyl.cal cyl_match_u cyl_match_v -s $thisrad >> $tree
	echo 0 2 1.5225225 1 bark_pat alias my_bark_mat bark_mat >> $tree
	echo my_bark_mat cone level$i 0 0 8 0 0 0 0 0 \
			$move $thisrad $lastrad >> $tree
	set spin=(`ev "rand($i)*360" "rand($i+$nlevels)*360" "rand($i+2*$nlevels)*360" "rand($i+3*$nlevels)*360" "rand($i+4*$nlevels)*360" "rand($i+5*$nlevels)*360" "rand($i+6*$nlevels)*360"`)
	xform -n b1 -s 1.1 -rz $spin[2] -ry -80 -rz $spin[1] -rz 5 -t 0 0 \
		`ev "$move*.42"` $oldtree >> $tree
	xform -n b2 -s 1.1 -rz $spin[3] -ry -78 -rz $spin[1] -rz 128 -t 0 0 \
		`ev "$move*.44"` $oldtree >> $tree
	xform -n b3 -s 1.1 -rz $spin[4] -ry -75 -rz $spin[1] -rz 255 -t 0 0 \
		`ev "$move*.40"` $oldtree >> $tree
	xform -n b4 -rz $spin[5] -ry -80 -rz $spin[1] -rz 58 -t 0 0 \
		`ev "$move*.92"` $oldtree >> $tree
	xform -n b5 -rz $spin[6] -ry -78 -rz $spin[1] -rz 181 -t 0 0 \
		`ev "$move*.84"` $oldtree >> $tree
	xform -n b6 -rz $spin[7] -ry -75 -rz $spin[1] -rz 297 -t 0 0 \
		`ev "$move*.88"` $oldtree >> $tree
	@ i++
end
#
# Send final tree
#
cat $tree
done:
rm -f $tree $oldtree
