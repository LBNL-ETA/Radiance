#!/bin/csh -f
# RCSid $Id: ra_pfm.csh,v 2.3 2004/12/30 18:33:50 greg Exp $
#
# Convert to/from Poskanzer Float Map image format using pvalue
#
if (`uname -p` == powerpc) then
	set machend=big
else
	set machend=little
endif
while ($#argv > 0)
	if ("$argv[1]" == "-r") then
		set reverse
	else if (! $?inp) then
		set inp="$argv[1]"
	else if (! $?out) then
		set out="$argv[1]"
	else
		goto userr
	endif
	shift argv
end
if ($?reverse) then
	if (! $?inp) then
		goto userr
	endif
	set opt=""
	set hl="`head -3 $inp:q`"
	if ("$hl[1]" == "Pf") then
		set opt=($opt -b)
	else if ("$hl[1]" != "PF") then
		echo "Input not a Poskanzer Float Map"
		exit 1
	endif
	if (`ev "if($hl[3],1,0)"`) then
		set filend=big
	else
		set filend=little
	endif
	if ($filend != $machend) then
		set opt=($opt -dF)
	else
		set opt=($opt -df)
	endif
	set res=($hl[2])
	if ($?out) then
		tail +4 $inp:q | pvalue -r -h $opt -y $res[2] +x $res[1] > $out:q
	else
		tail +4 $inp:q | pvalue -r -h $opt -y $res[2] +x $res[1]
	endif
	exit $status
endif
if (! $?inp) then
	goto userr
endif
set res=(`getinfo -d < $inp:q`)
if ($?out) then
	( echo PF ; echo $res[4] $res[2] ) > $out:q
	if ($machend == little) then
		echo "-1.000000" >> $out:q
	else
		echo "1.000000" >> $out:q
	endif
	pvalue -h -H -df $inp:q >> $out:q
else
	echo PF
	echo $res[4] $res[2]
	if ($machend == little) then
		echo "-1.000000"
	else
		echo "1.000000"
	endif
	pvalue -h -H -df $inp:q
endif
exit $status
userr:
echo "Usage: $0 input.pfm [output.hdr]"
echo "   or: $0 -r input.hdr [output.pfm]"
exit 1
