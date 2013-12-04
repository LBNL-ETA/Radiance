#!/usr/bin/perl
# RCSid $Id: objpict.pl,v 2.1 2013/12/04 20:07:07 greg Exp $
#
# Make a nice multi-view picture of an object
# Command line arguments contain materials and object files
#
# This is a re-write of Greg's objpict.csh and should be a drop-in
# replacement, with no funcionality added or removed

use strict;
use warnings;

use File::Temp qw/ tempdir /;
my $td = tempdir( CLEANUP => 1 );

my $xres = 1024;
my $yres = 1024;
my $rpict_cmd = "rpict -av .2 .2 .2 -x $xres -y $yres";

my $testroom = "$td/testroom.rad";
my $octree = "$td/op.oct";

# We need at least one Radiance file or a scene on STDIN (but not both)
if ($#ARGV < 0) {
	open(FH, ">$td/stdin.rad") or 
			die("objview: Can't write to temporary file $td/stdin.rad\n");
	while (<>) {
		print FH;
	}
	# Pretend stdin.rad was passed as argument.
	@ARGV = ("$td/stdin.rad");
}

# Create some lights and a box as back drop.
# The objects and view points will be inside the box.
open(FH, ">$testroom") or
		die("Can\'t write to temporary file $testroom");
print FH <<EndOfTestroom;
void plastic wall_mat  0  0  5  .681 .543 .686  0 .2
void light bright  0  0  3  3000 3000 3000

bright sphere lamp0  0  0  4  4 4 -4  .1
bright sphere lamp1  0  0  4  4 0 4  .1
bright sphere lamp2  0  0  4  0 4 4  .1

wall_mat polygon box.1540  0  0  12  5 -5 -5  5 -5 5  -5 -5 5  -5 -5 -5
wall_mat polygon box.4620  0  0  12  -5 -5 5  -5 5 5  -5 5 -5  -5 -5 -5
wall_mat polygon box.2310  0  0  12  -5 5 -5  5 5 -5  5 -5 -5  -5 -5 -5
wall_mat polygon box.3267  0  0  12  5 5 -5  -5 5 -5  -5 5 5  5 5 5
wall_mat polygon box.5137  0  0  12  5 -5 5  5 -5 -5  5 5 -5  5 5 5
wall_mat polygon box.6457  0  0  12  -5 5 5  -5 -5 5  5 -5 5  5 5 5
EndOfTestroom
close(FH);

my $dimstr = `getbbox -h @ARGV`;
chomp $dimstr;
# Values returned by getbbox are indented and delimited with multiple spaces.
$dimstr =~ s/^\s+//;   # remove leading spaces
my @dims = split(/\s+/, $dimstr);   # convert to array

# Find largest axes-aligned dimension
my @diffs = ($dims[1]-$dims[0], $dims[3]-$dims[2], $dims[5]-$dims[4]);
@diffs = reverse sort { $a <=> $b } @diffs;
my $size = $diffs[0];

# Define the four views
my $vw1 = "-vtl -vp 2 .5 .5 -vd -1 0 0 -vh 1 -vv 1";
my $vw2 = "-vtl -vp .5 2 .5 -vd 0 -1 0 -vh 1 -vv 1";
my $vw3 = "-vtl -vp .5 .5 2 -vd 0 0 -1 -vu -1 0 0 -vh 1 -vv 1";
my $vw4 = "-vp 3 3 3 -vd -1 -1 -1 -vh 20 -vv 20";

# Move objects so centre is at origin
my $xtrans = -1.0 * ($dims[0] + $dims[1]) / 2;
my $ytrans = -1.0 * ($dims[2] + $dims[3]) / 2;
my $ztrans = -1.0 * ($dims[4] + $dims[5]) / 2;
# Scale so that largest object dimension is unity
my $scale = 1 / $size;

my $cmd = "xform -t $xtrans $ytrans $ztrans -s $scale -t .5 .5 .5 @ARGV";
$cmd .= " |oconv $testroom - > $octree";
system "$cmd";

# Render four different views of the objects
system "$rpict_cmd $vw1 $octree > $td/right.hdr";
system "$rpict_cmd $vw2 $octree > $td/front.hdr";
system "$rpict_cmd $vw3 $octree > $td/down.hdr";
system "$rpict_cmd $vw4 $octree > $td/oblique.hdr";

# Compose the four views into one image
$cmd = "pcompos $td/down.hdr 0 $xres $td/oblique.hdr $xres $yres";
$cmd .= " $td/right.hdr 0 0 $td/front.hdr $xres 0";
$cmd .= " |pfilt -1 -r .6 -x /2 -y /2";
exec "$cmd";

#EOF
