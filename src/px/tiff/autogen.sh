#!/bin/sh
set -x
libtoolize --force --copy
aclocal -I .
autoheader
automake --foreign --add-missing --copy
autoconf

