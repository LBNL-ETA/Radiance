#!/bin/csh -f
# RCSid $Id: fieldcomb.csh,v 2.2 2005/09/02 04:24:37 greg Exp $
#
# Combine alternate lines in full frames for field rendering
#
# Expects numbered frames on command line, as given by ranimate
#
# If an odd number of frames is provided, the spare frame at the
# end is linked to $spare_name for the next run
#
# Written by Greg Ward for Iebele Atelier in August 2005
#
set spare_name=spare_fieldcomb_frame.pic
if ($#argv > 1) then
	if ("$argv[1]" == "-r") then
		set remove_orig
		shift argv
	endif
endif
if ($#argv < 2) then
	echo "Usage: $0 [-r] field1.pic field2.pic .."
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
	pcomb -e 'ro=ri(fld); go=gi(fld); bo=bi(fld)' \
		-e 'yd=yres-1-y; odd=.5*yd-floor(.5*yd)-.25' \
		-e 'fld=if(odd,1,2)' \
		$fields[$curfi]:q $fields[$nextfi]:q \
		> "${basenm}C$curfr.$ext"
	if ($?remove_orig) rm $fields[$curfi]:q $fields[$nextfi]:q
	@ curfr++
	@ curfi = $nextfi + 1
end
rm -f $spare_name
if ($curfi == $#fields) ln "${basenm}$curfi.$ext" $spare_name
