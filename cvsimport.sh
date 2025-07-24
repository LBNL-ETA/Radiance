#!/bin/bash
#
# This script imports the 'ray' module from the public Radiance CVS repository
# into a new Git repository.
#
# It requires an 'authors.txt' file in the same directory to map
# CVS usernames to Git author names and emails.

set -e

echo "Starting CVS import of Radiance 'ray' module..."

git cvsimport -k \
    -v \
    -o master \
    -d :pserver:anonymous@radiance-online.org:/home/cvsd/radiance \
    -A ../authors.txt \
    ray

echo "CVS import completed successfully."