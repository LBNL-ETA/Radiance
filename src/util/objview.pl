#!/usr/bin/perl
# RCSid $Id: objview.pl,v 2.3 2017/08/11 21:39:52 greg Exp $
#
# Make a nice view of an object
# Arguments are scene input files
#
# This is a re-write of Greg's original objview.csh.
# The only extra functionality is that we accept a scene on STDIN
# if no file name is given.
#
# Axel, Nov 2013

use strict;
use warnings;
use File::Temp qw/ tempdir /;

my $td = tempdir( CLEANUP => 1 );
my $octree = "$td/ov$$.oct";
my $lights = "$td/lt$$.rad";
my $rif = "$td/ov$$.rif";
my $ambf = "$td/af$$.amb";
my $raddev = "x11";   # default output device. Overwrite with -o
my $up = "Z";
my $vw = "XYZ";

my $opts = "";        # Options common to rad and glrad
my $rendopts = "";    # For render= line in rif file
my $usegl = 0;        # Run glrad instead of rad (Boolean).
my $radopt = 0;       # An option specific to rad was passed (Boolean).
my $glradopt = 0;     # An option specific to glrad was passed (Boolean).

while (@ARGV) {
	$_ = $ARGV[0];
	if (m/^-g\b/) {   # OpenGL output
		if ( $^O =~ /MSWin32/ ) {
			die("OpenGL view is not available under Windows.\n");
		}
		$usegl = 1;
	} elsif (m/^-u\b/) {   # up direction
		$up = $ARGV[1];
		shift @ARGV;
	} elsif ((m/^-s\b/) or (m/^-w/)) {   # silent, no warnings
		$opts .= " $_";
	} elsif (m/^-b/) {   # back face visibility
		$rendopts .= ' -bv';
	} elsif (m/^-v\b/) {   # standard view "[Xx]?[Yy]?[Zz]?[vlcahs]?"
		# Let rad do any error handling...
		$vw = $ARGV[1];
		shift @ARGV;
	} elsif (m/^-[nN]\b/) {   # No. of parallel processes
		$opts .= ' -N ' . $ARGV[1];
		$radopt = 1;
		shift @ARGV;
	} elsif (m/^-o\b/) {   # output device (rvu -devices)
		$raddev = $ARGV[1];
		$radopt = 1;
		shift @ARGV;
	} elsif ((m/^-V\b/) or (m/^-e\b/)) {   # print view, explicate variables
		# Think of those two as '-verbose'.
		$opts .= " $_";
		$radopt = 1;
	} elsif (m/^-S\b/) {   # full-screen stereo
		$opts .= " $_";
		$glradopt = 1;
	} elsif (m/^-\w/) {
		die("objview: Bad option: $_\n");
	} else {
		last;
	}
	shift @ARGV;
}

# We need at least one Radiance file or a scene on STDIN
if ($#ARGV < 0) {
	open(FH, ">$td/stdin.rad") or 
			die("objview: Can't write to temporary file $td/stdin.rad\n");
	while (<>) {
		print FH;
	}
	# Pretend stdin.rad was passed as argument.
	@ARGV = ("$td/stdin.rad");
}

# Make sure we don't confuse glrad and rad options.
if ($usegl) {
	if ($radopt) {
		die("objview: glrad output requested, but rad option passed.\n");
	}
} else {
	if ($glradopt) {
		die("objview: rad output requested, but glrad option passed.\n");
	}
}

open(FH, ">$lights") or
		die("objview: Can't write to temporary file $lights\n");
print FH <<EndOfLights;
void glow dim  0  0  4  .1 .1 .15  0
dim source background  0  0  4  0 0 1  360
void light bright  0  0  3  1000 1000 1000
bright source sun1  0  0  4  1 .2 1  5
bright source sun2  0  0  4  .3 1 1  5
bright source sun3  0  0  4  -1 -.7 1  5
EndOfLights
close(FH);

my @scenes = @ARGV;
push (@scenes, $lights);

# Make this work under Windoze
if ( $^O =~ /MSWin32/ ) {
	for my $i (0 .. $#scenes) {
		# rad doesn't like Windows-style backslashes.
		$scenes[$i] =~ s{\\}{/}g;
	}
	$octree =~ s{\\}{/}g;
	$ambf =~ s{\\}{/}g;
	$raddev = "qt";
}

my $scene = join(' ', @scenes);
open(FH, ">$rif") or
		die("objview: Can't write to temporary file $rif\n");
print FH <<EndOfRif;
scene= $scene
EXPOSURE= .5
UP= $up
view= $vw
OCTREE= $octree
oconv= -f
AMBF= $ambf
render= $rendopts
EndOfRif
close(FH);

if ($usegl) {
	system "glrad $opts $rif";
} else {
	system "rad -o $raddev $opts $rif";
}

#EOF
