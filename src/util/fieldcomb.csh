#!/bin/csh -f
# RCSid $Id: fieldcomb.csh,v 2.5 2005/09/05 16:13:08 greg Exp $
#
# Combine alternate lines in full frames for field rendering
#
# Expects numbered frames on command line, as given by ranimate
#
# If an odd number of frames is provided, the spare frame at the
# end is linked to $spare_name for the next run
#
# Written by Greg Ward for Iebele Abel in August 2005
#
set spare_name=spare_fieldcomb_frame.pic
set odd_first=0
while ($#argv > 1)
	switch ($argv[1])
	case -r*:
		set remove_orig
		breaksw
	case -o*:
		set odd_first=1
		breaksw
	case -e*:
		set odd_first=0
		breaksw
	case -f*:
		shift argv
		set outfile="$argv[1]"
		breaksw
	default:
		if ("$argv[1]" !~ -*) break
		echo "Unknown option: $argv[1]"
		exit 1
	endsw
	shift argv
end
if ($#argv < 2 || ($?outfile && $#argv > 2)) then
	echo "Usage: $0 [-e|-o][-r] [-f combined.pic] field1.pic field2.pic .."
	exit 1
endif
set f1=$argv[1]:q
set ext=$f1:e
set basenm="`echo $f1:q | sed 's/[0-9]*\.'$ext'//'`"
set curfi=`echo $f1:q | sed 's/^[^1-9]*\(.[0-9]*\)\.'$ext'$/\1/'`
set fields=($argv[*]:q)
if (-r $spare_name) then
	set fields=($spare_name $fields:q)
	@ curfi--
endif
@ curfr = $curfi / 2
set curfi=1
while ($curfi < $#fields)
	@ nextfi = $curfi + 1
	if ($curfr < 10) then
		set fid=000$curfr
	else if ($curfr < 100) then
		set fid=00$curfr
	else if ($curfr < 1000) then
		set fid=0$curfr
	else
		set fid=$curfr
	endif
       	set outf="${basenm}C$fid.$ext"
	if ($?outfile) set outf=$outfile:q
	pcomb -e 'ro=ri(fld); go=gi(fld); bo=bi(fld)' \
		-e 'yd=yres-1-y; odd=.5*yd-floor(.5*yd)-.25' \
		-e "fld=if(odd,2-$odd_first,1+$odd_first)" \
		$fields[$curfi]:q $fields[$nextfi]:q \
		> $outf:q
	if ($?remove_orig) rm $fields[$curfi]:q $fields[$nextfi]:q
	@ curfr++
	@ curfi = $nextfi + 1
end
rm -f $spare_name
if ($curfi == $#fields) ln "${basenm}$curfi.$ext" $spare_name
