#!/usr/bin/perl
# RCSid $Id: bsdfview.pl,v 2.4 2017/08/26 16:07:22 greg Exp $
#
# Call bsdf2rad to render BSDF and start viewing it.
# Arguments are BSDF XML or SIR file(s)
#
use strict;
use warnings;
use File::Temp qw/ tempdir /;

my $td = tempdir( CLEANUP => 0 );
my $octree = "$td/bv$$.oct";
my $ambf = "$td/af$$.amb";
my $raddev = "x11";	# default output device. Overwrite with -o
my $qual = "Low";
my $usetrad = 0;

my $opts = "";	# Options common to rad
my $rendopts = "-w-";	# For render= line in rif file

while (@ARGV) {
	$_ = $ARGV[0];
	if ((m/^-s\b/) or (m/^-w/)) {	# silent, no warnings
		$opts .= " $_";
	} elsif (m/^-v\b/) {	# standard view
		# Let rad do any error handling...
		$opts .= qq( -v "$ARGV[1]");
		shift @ARGV;
	} elsif (m/^-[nN]\b/) {	# No. of parallel processes
		$opts .= ' -N ' . $ARGV[1];
		shift @ARGV;
	} elsif (m/^-o\b/) {	# output device (rvu -devices)
		$raddev = $ARGV[1];
		shift @ARGV;
	} elsif (m/^-q/) {	# quality setting
		$qual = $ARGV[1];
		shift @ARGV;
	} elsif ((m/^-V\b/) or (m/^-e\b/)) {   # print view, explicate variables
		# Think of those two as '-verbose'.
		$opts .= " $_";
	} elsif (m/^-t\b/) {	# start trad instead of rad
		$usetrad = 1;
	} elsif (m/^-\w/) {
		die("bsdfview: Bad option: $_\n");
	} else {
		last;
	}
	shift @ARGV;
}

# We need at least one XML or SIR file
if ($#ARGV < 0) {
	die("bsdfview: missing input XML or SIR file(s)\n");
}

if (length($opts) and $usetrad) {
	die("bsdfview: rad options not supported when calling trad (-t)\n");
}

my @objects = @ARGV;

# Make this work under Windoze
if ( $^O =~ /MSWin32/ ) {
	for my $i (0 .. $#objects) {
		# rad doesn't like Windows-style backslashes.
		$objects[$i] =~ s{\\}{/}g;
	}
	$octree =~ s{\\}{/}g;
	$ambf =~ s{\\}{/}g;
	$raddev = "qt";
}

my $name = $objects[0];
$name =~ s{.*/}{};		# remove leading path
$name =~ s{\.[^.]+$}{};		# remove file extension

my $rif = "$name.rif";

if (-e $rif) {			# RIF already exists?
	print "Attempting to run with existing rad input file '$rif'\n";
	if ($usetrad) {
		system "trad $rif";
	} else {
		system "rad -o $raddev -w $opts $rif QUA=$qual";
	}
	die("\nTry removing '$rif' and starting again\n\n") if $?;
	exit;
}

print "bsdfview: creating rad input file '$rif'\n";

my $scene = qq("!bsdf2rad @objects");		# let bsdf2rad do complaining

my $objects = join(' ', @objects);
open(FH, ">$rif") or
		die("bsdfview: Can't write to temporary file $rif\n");
print FH <<EndOfRif;
scene= $scene
objects= $objects
ZONE= E -35 35 -20 15 -5 15
PICTURE= $name
RESOLU= 1024
EXPOSURE= 1
UP= +Z
OCTREE= $octree
oconv= -w -f
AMBF= $ambf
QUAL= $qual
render= $rendopts
view= def -vp 0 -50 50 -vd 0 50 -50 -vh 45 -vv 30
view= fr -vp 15 -30 30 -vd 0 30 -30
view= br -vp -15 -30 30 -vd 0 30 -30
view= ft -vta -vp 15 0 0 -vd 0 0 1 -vu 0 1 0 -vh 200 -vv 200
view= bt -vta -vp -15 0 0 -vd 0 0 1 -vu 0 1 0 -vh 200 -vv 200
view= pr -vtl -vp 0 0 20 -vd 0 0 -1 -vu 0 1 0 -vv 35 -vh 65
view= pt -vtl -vp 0 0 -10 -vd 0 0 1 -vu 0 1 0 -vv 35 -vh 65
EndOfRif
close(FH);

if ($usetrad) {
	system "trad $rif";
} else {
	system "rad -o $raddev $opts $rif";
}
#EOF
