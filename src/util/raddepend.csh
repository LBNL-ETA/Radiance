#!/bin/csh -f
# RCSid: $Id: raddepend.csh,v 2.7 2003/02/22 02:07:30 greg Exp $
#
# Find scene dependencies in this directory
#
set es=1
onintr quit
rm -f EMPTY
echo -n > EMPTY
sleep 2
( ls $* | sed -e 's~/~\\/~g' -e 's@^@/^@' -e 's@$@$/d@' ; echo '/^EMPTY$/,$d' ) > /tmp/sed$$
getbbox -w $* >/dev/null
set es=$status
if ( $es == 0 ) then
	sync
	sleep 2
	ls -tuL | sed -f /tmp/sed$$ | sort
endif
quit:
rm -f /tmp/sed$$ EMPTY
exit $es
