#!/bin/csh -f
# RCSid $Id: ra_pfm.csh,v 2.1 2004/11/25 14:47:03 greg Exp $
#
# Convert to/from Poskanzer Float Map image format using pvalue
#
if (`uname -p` == powerpc) then
	set format="-dF"
else
	set format="-df"
endif
while ($#argv > 0)
	if ("$argv[1]" == "-r") then
		set reverse
	else
		set inp="$argv[1]"
	endif
	shift argv
end
if ($?reverse) then
	if (! $?inp) then
		goto userr
	endif
	set hl="`head -3 $inp:q`"
	if ("$hl[1]" != "PF") then
		echo "Input not a Poskanzer Float Map"
		exit 1
	endif
	set res=($hl[2])
	tail +4 $inp:q | pvalue -r -h -y $res[1] +x $res[2] $format
	exit $status
endif
if (! $?inp) then
	goto userr
endif
set res=(`getinfo -d < $inp:q`)
echo PF
echo $res[4] $res[2]
echo "-1.000000"
pvalue -h -H $format $inp:q
exit $status
userr:
echo "Usage: $0 input.pfm > output.hdr"
echo "   or: $0 -r input.hdr > output.pfm"
exit 1
