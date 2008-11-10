#!/bin/csh -f
# RCSid: $Id: ivpict.csh,v 2.3 2008/11/10 19:08:19 greg Exp $
#
# Render requested Instant View
#
if ( $#argv < 1 ) then
	echo Usage: $0 directory '[options]'
	exit 1
endif
set dir=$argv[1]
set opt="$argv[2-]"
exec pinterp $opt $dir/pvw_{x,X,y,Y,z,Z}l.{hdr,z}
