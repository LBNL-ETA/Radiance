#!/bin/csh -f
#
# Apply Reinhard's simple tone operator
#
if ($#argv != 3) then
	echo Usage: $0 key input.hdr output.tif
	exit 1
endif
set lavg=`pvalue -h -H -b -d $argv[2]:q | total -p -m`
set lmax=`pextrem $argv[2]:q | sed -n '2s/^[1-9][0-9]* [1-9][0-9]*	[^ ]* \([^ ]*\) .*$/\1/p'`
pcomb -e "Lavg:$lavg;Lwht:$lmax;a:$argv[1]" -f reinhard.cal $argv[2]:q \
	| ra_tiff - $argv[3]:q
