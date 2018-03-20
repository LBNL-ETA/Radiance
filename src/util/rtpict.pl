#!/usr/bin/perl -w
# RCSid $Id: rtpict.pl,v 2.4 2018/03/20 23:05:48 greg Exp $
#
# Run rtrace in parallel mode to simulate rpict -n option
#
#	G. Ward
#
use strict;
# we'll call rpict if this is not overridden
my $nprocs = 1;
# rtrace options and the associated number of arguments
my %rtraceC = ("-dt",1, "-dc",1, "-dj",1, "-ds",1, "-dr",1, "-dp",1,
		"-ss",1, "-st",1, "-e",1, "-am",1,
		"-ab",1, "-af",1, "-ai",1, "-aI",1, "-ae",1, "-aE",1,
		"-av",3, "-aw",1, "-aa",1, "-ar",1, "-ad",1, "-as",1,
		"-me",3, "-ma",3, "-mg",1, "-ms",1, "-lr",1, "-lw",1);
# boolean rtrace options
my @boolO = ("-w", "-bv", "-dv", "-i", "-u");
# view options and the associated number of arguments
my %vwraysC = ("-vf",1, "-vtv",0, "-vtl",0, "-vth",0, "-vta",0, "-vts",0, "-vtc",0,
		"-x",1, "-y",1, "-vp",3, "-vd",3, "-vu",3, "-vh",1, "-vv",1,
		"-vo",1, "-va",1, "-vs",1, "-vl",1, "-pa",1, "-pj",1, "-pd",1);
# options we need to silently ignore
my %ignoreC = ("-t",1, "-ps",1, "-pt",1, "-pm",1);
# Starting options for rtrace (rpict values)
my @rtraceA = split(' ', "rtrace -ffc -u- -dt .05 -dc .5 -ds .25 -dr 1 -aa .2 -ar 64 -ad 512 -as 128 -lr 7 -lw 1e-03");
my @vwraysA = ("vwrays", "-ff", "-pj", ".67");
my @vwrightA = ("vwright", "-vtv");
my @rpictA = ("rpict");
my $outpic;
OPTION:				# sort through options
while ($#ARGV >= 0 && "$ARGV[0]" =~ /^[-\@]/) {
	# Check for file inclusion
	if ("$ARGV[0]" =~ /^\@/) {
		open my $handle, '<', substr($ARGV[0], 1);
		shift @ARGV;
		chomp(my @args = <$handle>);
		close $handle;
		unshift @ARGV, split(/\s+/, "@args");
		next OPTION;
	}
	# Check booleans
	for my $boopt (@boolO) {
		if ("$ARGV[0]" =~ ('^' . $boopt . '[-+01tfynTFYN]$')) {
			push @rtraceA, $ARGV[0];
			push @rpictA, shift(@ARGV);
			next OPTION;
		}
	}
	# Check view options
	if (defined $vwraysC{$ARGV[0]}) {
		push @vwraysA, $ARGV[0];
		my $isvopt = ("$ARGV[0]" =~ /^-v/);
		push(@vwrightA, $ARGV[0]) if ($isvopt);
		push @rpictA, shift(@ARGV);
		for (my $i = $vwraysC{$rpictA[-1]}; $i-- > 0; ) {
			push @vwraysA, $ARGV[0];
			push(@vwrightA, $ARGV[0]) if ($isvopt);
			push @rpictA, shift(@ARGV);
		}
		next OPTION;
	}
	# Check rtrace options
	if (defined $rtraceC{$ARGV[0]}) {
		push @rtraceA, $ARGV[0];
		push @rpictA, shift(@ARGV);
		for (my $i = $rtraceC{$rpictA[-1]}; $i-- > 0; ) {
			push @rtraceA, $ARGV[0];
			push @rpictA, shift(@ARGV);
		}
		next OPTION;
	}
	# Check options to ignore
	if (defined $ignoreC{$ARGV[0]}) {
		push @rpictA, shift(@ARGV);
		for (my $i = $ignoreC{$rpictA[-1]}; $i-- > 0; ) {
			push @rpictA, shift(@ARGV);
		}
		next OPTION;
	}
	# Check for output file or number of processes
	if ("$ARGV[0]" eq "-o") {
		shift @ARGV;
		$outpic = shift(@ARGV);
	} elsif ("$ARGV[0]" eq "-n") {
		shift @ARGV;
		$nprocs = shift(@ARGV);
	} else {
		die "Unsupported option: " . $ARGV[0] . "\n";
	}
}
die "Number of processes must be positive" if ($nprocs <= 0);
if (defined $outpic) {		# redirect output?
	die "File '$outpic' already exists\n" if (-e $outpic);
	open STDOUT, '>', "$outpic";
}
if ($nprocs == 1) {		# may as well run rpict?
	push(@rpictA, $ARGV[0]) if ($#ARGV == 0);
	exec @rpictA ;
}
die "Need single octree argument\n" if ($#ARGV != 0);
push @rtraceA, (`@vwraysA -d`);
chomp $rtraceA[-1];
push @rtraceA, ("-n", "$nprocs");
push @rtraceA, $ARGV[0];
my @view = (`@vwrightA 0`);
exec qq{@vwraysA | @rtraceA | getinfo -a "VIEW=@view"};
