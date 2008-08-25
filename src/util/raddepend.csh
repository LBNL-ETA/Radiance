#!/bin/csh -f
# RCSid: $Id: raddepend.csh,v 2.8 2008/08/25 04:50:32 greg Exp $
#
# Find scene dependencies in this directory
#
set es=1
onintr quit
rm -f EMPTY
echo -n > EMPTY
sleep 2
set sedf=`mktemp /tmp/sed.XXXXXX`
( ls $* | sed -e 's~/~\\/~g' -e 's@^@/^@' -e 's@$@$/d@' ; echo '/^EMPTY$/,$d' ) > $sedf
getbbox -w $* >/dev/null
set es=$status
if ( $es == 0 ) then
	sync
	sleep 2
	ls -tuL | sed -f $sedf | sort
endif
quit:
rm -f $sedf EMPTY
exit $es
