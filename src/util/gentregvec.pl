#!/usr/bin/perl -w
# RCSid $Id: gentregvec.pl,v 2.2 2009/06/16 17:18:38 greg Exp $
#
# Generate Tregenza vector for a given sky description
#
#	G. Ward
#
use strict;
$ENV{RAYPATH} = ':/usr/local/lib/ray' if (! $ENV{RAYPATH});
my $tregsamp = "tregsamp.dat";
my @skycolor = (0.960, 1.004, 1.118);
while ($#ARGV >= 0) {
	if ("$ARGV[0]" eq "-c") {
		@skycolor = @ARGV[1..3];
		shift @ARGV; shift @ARGV; shift @ARGV;
	}
	shift @ARGV;
}
# Load sky description into line array, separating sun if one
my @skydesc;
my $lightline;
my @sunval;
my $sunline;
my $srcmod;	# putting this inside loop breaks code(?!)
while (<>) {
	push @skydesc, $_;
	if (/^\w+\s+light\s+/) {
		s/\s+$//; s/^.*\s//;
		$srcmod = $_;
		$lightline = $#skydesc;
	} elsif (defined($srcmod) && /^($srcmod)\s+source\s/) {
		@sunval = split(/\s+/, $skydesc[$lightline + 3]);
		shift @sunval;
		$sunline = $#skydesc;
	}
}
die "Empty input!\n" if (! @skydesc);
# Strip out the solar source if present
my @sundir;
if (defined $sunline) {
	@sundir = split(/\s+/, $skydesc[$sunline + 3]);
	shift @sundir;
	undef @sundir if ($sundir[2] <= 0);
	splice(@skydesc, $sunline, 5);
}
# Get Tregenza sample file
my $found;
foreach my $dir (split /:/, $ENV{RAYPATH}) {
	if (-r "$dir/$tregsamp") {
		$found = "$dir/$tregsamp";
		last;
	}
}
if (! $found) {
	print STDERR "Cannot find file '$tregsamp'\n";
	exit 1;
}
$tregsamp = $found;
# Create octree for rtrace
my $octree = "/tmp/gtv$$.oct";
open OCONV, "| oconv - > $octree";
print OCONV @skydesc;
print OCONV "skyfunc glow skyglow 0 0 4 @skycolor 0\n";
print OCONV "skyglow source sky 0 0 4 0 0 1 360\n";
close OCONV;
# Run rtrace and average output for every 64 samples
my @tregval = `rtrace -h -faf -ab 0 -w < $tregsamp $octree | total -if3 -m -64`;
if ($#tregval != 145) {
	print STDERR "Expected 9344 samples in $tregsamp -- found ",
			`wc -l < $tregsamp`, "\n";
	exit 1;
}
unlink $octree;
# Find closest 3 patches to sun and divvy up direct solar contribution
my @bestdir;
if (@sundir) {
	my $somega = ($sundir[3]/360)**2 * 3.141592654**3;
	my @tindex = (30, 60, 84, 108, 126, 138, 144, 145);
	my @tomega = (0.0435449227, 0.0416418006, 0.0473984151,
			0.0406730411, 0.0428934136, 0.0445221864,
			0.0455168385, 0.0344199465);
	my $cmd = "total -m -64 < $tregsamp | tail +2 | rcalc " .
		q{-e 'Dx=$4;Dy=$5;Dz=$6' } .
		"-e 'dot=Dx*$sundir[0] + Dy*$sundir[1] + Dz*$sundir[2]' " .
		q{-e '$1=acos(dot);$2=recno' | sort -n};
	@bestdir = `$cmd | head -3`;
	my @ang; my @ndx;
	($ang[0], $ndx[0], $ang[1], $ndx[1], $ang[2], $ndx[2]) =
		(	split(/\s+/, $bestdir[0]),
			split(/\s+/, $bestdir[1]),
			split(/\s+/, $bestdir[2]) );
	my $wtot = 1/($ang[0]+.02) + 1/($ang[1]+.02) + 1/($ang[2]+.02);
	for my $i (0..2) {
		my $ti = 0;
		while ($ndx[$i] > $tindex[$ti]) { $ti++ }
		my $wt = $wtot/($ang[$i]+.02) * $somega / $tomega[$ti];
		my @scolor = split(/\s+/, $tregval[$ndx[$i]]);
		for my $j (0..2) { $scolor[$j] += $wt * $sunval[$j]; }
		$tregval[$ndx[$i]] = "@scolor\n";
	}
}
# Output our final vector
print @tregval;
