#!/bin/sh
# RCSid: $Id: getviews.sh,v 1.2 2021/02/19 16:04:44 greg Exp $
sed -n 's/^view *= *\([a-zA-Z][a-zA-Z0-9]*\) .*$/\1/p' $*
