#!/bin/sh
sed -n 's/^view *= *\([a-zA-Z][a-zA-Z0-9]*\) .*$/\1/p' $*
