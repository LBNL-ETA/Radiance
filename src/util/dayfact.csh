#!/bin/csh -f
# RCSid: $Id: dayfact.csh,v 2.5 2003/02/22 02:07:30 greg Exp $
#
# Interactive script to calculate daylight factors
#
set nofile="none"
set illumpic=$nofile
set octree=$nofile
set dfpict=$nofile
set ilpict=$nofile
set dspict=$nofile
set nodaysav=1
set fcopts=($*)
set designlvl=500
set wporig=(0 0 0)
set wpsize=(1 1)
set rtargs=(-ab 1 -ad 256 -as 128 -aa .15 -av .3 .3 .3)
set maxres=128

alias readvar 'echo -n Enter \!:1 "[$\!:1]: ";set ans="$<";if("$ans" != "")set \!:1="$ans"'

cat <<_EOF_
			DAYLIGHT FACTOR CALCULATION

This script calculates daylight factors, illuminance levels, and/or
energy savings due to daylight on a rectangular workplane and produces
a contour plot from the result.  The input is a Radiance scene description
(and octree) and the output is one or more color Radiance pictures.

_EOF_
echo "Have you already calculated an illuminance picture using dayfact?"
readvar illumpic
if ( $illumpic != $nofile ) then
	if ( ! -r $illumpic ) then
		echo "Cannot read $illumpic"
		exit 1
	endif
	set title=$illumpic:r
	set gotillumpic
	goto getgenskyf
endif

readvar octree
if ( $octree == $nofile || ! -f $octree ) then
	echo "You must first create an octree with"
	echo "oconv before running this script."
	exit 1
endif
set title="$octree:r"
echo -n "Is the z-axis your zenith direction? "
if ( "$<" !~ [yY]* ) then
	echo "I'm sorry, you cannot use this script"
	exit 1
endif
echo "What is the origin (smallest x y z coordinates) of the workplane?"
readvar wporig
set wporig=($wporig)
echo "What is the x and y size (width and length) of the workplane?"
readvar wpsize
set wpsize=($wpsize)
set wpres=(`rcalc -n -e '$1=if(l,'"floor($maxres*$wpsize[1]/$wpsize[2]),$maxres);"'$2=if(l,'"$maxres,floor($maxres*$wpsize[2]/$wpsize[1]));l=$wpsize[2]-$wpsize[1]"`)
set rtargs=($rtargs -ar `getinfo -d<$octree|rcalc -e '$1=floor(16*$4/'"($wpsize[1]+$wpsize[2]))"`)
echo "What calculation options do you want to give to rtrace?"
echo "(It is very important to set the -a* options correctly.)"
readvar rtargs
echo "Do you want to save the illuminance picture for later runs?"
readvar illumpic
############
getgenskyf:
set genskyf=$nofile
echo "In what scene file is the gensky command located?"
readvar genskyf
if ( $genskyf == $nofile ) then
	echo "You will not be able to compute daylight factors"
	echo "or energy savings since there is no gensky file."
else
	xform -e $genskyf > /usr/tmp/gsf$$
	grep '^# gensky ' /usr/tmp/gsf$$
	if ( $status ) then
		echo "The file $genskyf does not contain a gensky command\!"
		rm -f /usr/tmp/gsf$$
		goto getgenskyf
	endif
	set title=$title\ `sed -n 's/^# gensky  *\([0-9][0-9]*  *[0-9][0-9]*  *[0-9][0-9.]*\).*$/\1/p' /usr/tmp/gsf$$`
	set extamb=`sed -n 's/^# Ground ambient level: //p' /usr/tmp/gsf$$`
	grep -s '^# gensky .* -c' /usr/tmp/gsf$$
	set nodaysav=$status
	rm -f /usr/tmp/gsf$$
	if ( $nodaysav ) then
		echo "The gensky command was not done for an overcast sky"
		echo "(-c option), so energy savings cannot be calculated."
		echo -n "Continue anyway? "
		if ( "$<" =~ [nN]* ) then
			exit 0
		endif
	endif
endif
echo "Illuminance contour picture if you want one"
readvar ilpict
if ( $?extamb ) then
	echo "Daylight factor contour picture if you want one"
	readvar dfpict
endif
if ( ! $nodaysav ) then
	echo "Energy savings contour picture if you want one"
	readvar dspict
	if ( $dspict != $nofile ) then
		echo "Workplane design level (lux)"
		readvar designlvl
	endif
endif
if ( $ilpict == $nofile && $dfpict == $nofile && $dspict == $nofile ) then
	echo "Since you don't want any output, I guess we're done."
	exit 0
endif
echo "Title for output picture"
readvar title
set sctemp=/usr/tmp/sc$$.csh
cat <<'_EOF_' > $sctemp
if ( $illumpic != $nofile ) then
	set iltemp=""
else
	set iltemp=/usr/tmp/il$$.pic
	set illumpic=$iltemp
endif
set tltemp=/usr/tmp/tl$$.pic
set dstemp=/usr/tmp/ds$$.pic
set temp1=/usr/tmp/tfa$$
set tempfiles=($iltemp $sctemp $tltemp $dstemp $temp1)
echo "Your dayfact job is finished."
echo "Please check for error messages below."
echo ""
set echo
if ( ! $?gotillumpic ) then
	cnt $wpres[2] $wpres[1] \
	| rcalc -e '$1=($2+.5)/'"$wpres[1]*$wpsize[1]+$wporig[1]" \
		-e '$2=(1-($1+.5)/'"$wpres[2])*$wpsize[2]+$wporig[2]" \
		-e '$3='"$wporig[3]" -e '$4=0;$5=0;$6=1' \
	| rtrace $rtargs -h+ -I+ -ov -fac -x $wpres[1] -y $wpres[2] $octree \
	> $temp1
	pfilt -h 20 -n 0 -x 300 -y 300 -p 1 -r 1 $temp1 > $illumpic
endif
set maxval=`getinfo < $illumpic | rcalc -i 'EXPOSURE=${e}' -e '$1=3/e'`
psign -h 42 " $title " | pfilt -1 -x /2 -y /2 > $tltemp
'_EOF_'
if ( $ilpict != $nofile ) then
	echo 'falsecolor -cb -l Lux -s "$maxval*179" \\
		$fcopts -m 179 -ip $illumpic \\
		| pcompos -a 1 - $tltemp > $ilpict' >> $sctemp
endif
if ( $dfpict != $nofile ) then
	echo 'falsecolor -cb -l DF -s 50 \\
		$fcopts -m "100/PI/$extamb" -ip $illumpic \\
		| pcompos -a 1 - $tltemp > $dfpict' >> $sctemp
endif
if ( $dspict != $nofile ) then
	echo 'pcomb -e "lo=1-$designlvl/20000*3.1416*$extamb/li(1)" \\
		-o $illumpic | falsecolor -cb -l "%Save" -s 100 \\
		$fcopts -m 100 -p $illumpic \\
		| pcompos -a 1 - $tltemp > $dspict' >> $sctemp
endif
echo 'rm -f $tempfiles' >> $sctemp
(source $sctemp) |& mail `whoami` &
echo "Your job is started in the background."
echo "You will be notified by mail when it is done."
