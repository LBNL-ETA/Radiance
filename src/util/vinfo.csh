#!/bin/csh -f
# RCSid: $Id: vinfo.csh,v 3.2 2004/11/19 23:00:49 greg Exp $
#
# Edit information header in Radiance file
#
set plist=()
set ilist=()
foreach f ($*)
	if (! -f $f:q) then
		echo "${f}: no such file or directory"
		continue
	endif
	if (! -w $f:q) then
		echo "$f is read-only"
		sleep 1
		continue
	endif
	set info="$f.info"
	getinfo < $f:q > $info:q
	set plist=($plist:q $f:q)
	set ilist=($ilist:q $info:q)
end
vi $ilist:q
set i=1
while ( $i <= $#plist )
	set f=$plist[$i]:q
	set info=$ilist[$i]:q
	if ("`tail -1 $info:q`" != "") then
		echo "" >> $info:q
	endif
	getinfo < $f:q | cmp -s - $info:q
	if ($status != 0) then
		getinfo - < $f:q >> $info:q
		mv $info:q $f:q
	else
		rm $info:q
	endif
	@ i++
end
