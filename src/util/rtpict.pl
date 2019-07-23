#!/usr/bin/perl -w
# RCSid $Id: rtpict.pl,v 2.8 2019/07/23 01:03:30 greg Exp $
#
# Run rtrace in parallel mode to simulate rpict -n option
# May also be used to render layered images with -o* option
#
#	G. Ward
#
use strict;
# we'll run rpict if no -n or -o* option
my $nprocs = 1;
# rtrace options and the associated number of arguments
my %rtraceC = ("-dt",1, "-dc",1, "-dj",1, "-ds",1, "-dr",1, "-dp",1,
		"-ss",1, "-st",1, "-e",1, "-am",1,
		"-ab",1, "-af",1, "-ai",1, "-aI",1, "-ae",1, "-aE",1,
		"-av",3, "-aw",1, "-aa",1, "-ar",1, "-ad",1, "-as",1,
		"-me",3, "-ma",3, "-mg",1, "-ms",1, "-lr",1, "-lw",1);
# boolean rtrace options
my @boolO = ("-w", "-bv", "-dv", "-i", "-u");
my $irrad = 0;
# view options and the associated number of arguments
my %vwraysC = ("-vf",1, "-vtv",0, "-vtl",0, "-vth",0, "-vta",0, "-vts",0, "-vtc",0,
		"-x",1, "-y",1, "-vp",3, "-vd",3, "-vu",3, "-vh",1, "-vv",1,
		"-vo",1, "-va",1, "-vs",1, "-vl",1, "-pa",1, "-pj",1, "-pd",1);
# options we need to silently ignore
my %ignoreC = ("-t",1, "-ps",1, "-pt",1, "-pm",1);
# Starting options for rtrace (rpict values)
my @rtraceA = split(' ', "rtrace -u- -dt .05 -dc .5 -ds .25 -dr 1 " .
				"-aa .2 -ar 64 -ad 512 -as 128 -lr 7 -lw 1e-03");
my @vwraysA = ("vwrays", "-ff", "-pj", ".67");
my @vwrightA = ("vwright", "-vtv");
my @rpictA = ("rpict");
my @rsplitA = ("rsplit", "'-t	'", "-ih", "-iH", "-oh", "-oH");
my $refDepth="";
# Supported rtrace -o* options and data type encoders
my %rtoutC = (
	v =>	"radiance.hdr",
	r =>	"mirrored.hdr",
	x =>	"unmirrored.hdr",
	i =>	"irradiance.hdr",
	l =>	"effective.dpt",
	L =>	"firstsurface.dpt",
	R =>	"mirrored.dpt",
	X =>	"unmirrored.dpt",
	n =>	"perturbed.nrm",
	N =>	"unperturbed.nrm",
	s =>	"surface.idx",
	m =>	"modifier.idx",
	M =>	"material.idx"
);
my $outpatt = "^-o[vrxlLRXnNsmM]+";
# Supported compressed format suffixes and compressors, decompressors
my %comprC = (
	".tgz",	["tar czf ", "tar xzf "],
	".zip",	["zip -q -r ", "unzip -q "]
);
my $outlyr;
my $outdir;
my $outpic;
my $outzbf;
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
		if ("$ARGV[0]" =~ ('^' . $boopt . '[-+01tfynTFYN]?$')) {
			push @rtraceA, $ARGV[0];
			push @rpictA, shift(@ARGV);
			if ("$ARGV[0]" eq "-i") {
				$irrad = ! $irrad;
			} elsif ("$ARGV[0]" =~ /^-i[-+01tfynTFYN]/) {
				$irrad = ("$ARGV[0]" =~ /^-i[+1tyTY]/);
			}
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
	} elsif ("$ARGV[0]" eq "-z") {
		push @rpictA, shift(@ARGV);
		$outzbf = $ARGV[0];
		push @rpictA, shift(@ARGV);
	} elsif ("$ARGV[0]" eq "-n") {
		shift @ARGV;
		$nprocs = shift(@ARGV);
	} elsif ("$ARGV[0]" eq "-d") {
		shift @ARGV;
		$refDepth = " -d " . shift(@ARGV);
	} elsif ("$ARGV[0]" =~ "$outpatt") {
		push @rtraceA, $ARGV[0];
		$outlyr = substr(shift(@ARGV), 2);
		$outdir = shift(@ARGV);
	} else {
		die "Unsupported option: " . $ARGV[0] . "\n";
	}
}
my %rcodeC = (
	".hdr",	["-of3", "!pvalue -r -df -u"],
	".dpt",	["-of", "!rcode_depth$refDepth -ff"],
	".nrm",	["-of3", "!rcode_norm -ff"],
	".idx",	["-oa", '!rcode_ident "-t	"']
);
die "Number of processes must be positive" if ($nprocs <= 0);
if (defined $outdir) {		# check conflicting options
	die "Options -o and -o* are mutually exclusive\n" if (defined $outpic);
	die "Options -z and -o* are mutually exclusive\n" if (defined $outzbf);
}
if (defined $outpic) {		# redirect output?
	die "File '$outpic' already exists\n" if (-e $outpic);
	open STDOUT, '>', "$outpic";
}
				# may as well run rpict?
if ($nprocs == 1 && !defined($outdir)) {
	push(@rpictA, $ARGV[0]) if ($#ARGV == 0);
	exec @rpictA ;
}
push @rtraceA, ("-n", "$nprocs");
die "Need single octree argument\n" if ($#ARGV != 0);
my $oct = $ARGV[0];
my $view = `@vwrightA 0`;
my @res = split(/\s/, `@vwraysA -d`);
if (defined $outdir) {		# layered output directory?
	my $pkg;
	my ($pkgext) = ($outdir =~ /(\.[^.]+)$/);
	
	if ($pkgext) {
		die "Unknown compressor for '$pkgext'\n" if (!defined $comprC{$pkgext});
		$pkg = $outdir;
		die "Output '$pkg' exists as directory!\n" if (-d $pkg);
		$outdir = substr($pkg, 0, length($pkg)-length($pkgext));
		if (-e $outdir) {
			die "Both '$outdir' and '$pkg' exist!\n" if (-e $pkg);
		} elsif (-e $pkg) {
			system "$comprC{$pkgext}[1] '$pkg'";
			die "Cannot decompress '$pkg'\n" if ( $? );
		}
	}
	if (! -d $outdir) {
		mkdir $outdir or die("Cannot create output directory '$outdir'\n");
	}
				# construct rsplit command
	foreach (0..(length($outlyr) - 1)) {
		my $oval=substr($outlyr, $_, 1);
		$oval = "i" if ($oval eq "v" && $irrad);
		die "Duplicate letter in -o$outlyr\n" if (!defined $rtoutC{$oval});
		my ($otyp) = ($rtoutC{$oval} =~ /(\.[^.]+)$/);
		push @rsplitA, $rcodeC{$otyp}[0];
		push @rsplitA, qq{'$rcodeC{$otyp}[1] > "$outdir/$rtoutC{$oval}"'};
		delete $rtoutC{$oval};
	}
	system "@vwraysA | @rtraceA -fff @res $oct | " .
			"getinfo -a 'VIEW=$view' | @rsplitA";
	die "Error running rtrace\n" if ( $? );
	if (defined $pkg) {	# compress folder if requested
		unlink $pkg;
		system "$comprC{$pkgext}[0] '$pkg' '$outdir'";
		die "Cannot compress folder '$outdir'\n" if ( $? );
		system "rm -r '$outdir'";
	}
	exit;			# all done with layered output
}
if (defined $outzbf) {		# generating depth buffer?
	exec "@vwraysA | @rtraceA -fff -olv @res $oct | " .
		"rsplit -ih -iH -f " .
			"-of $outzbf " .
			"-oh -oH -of3 '!pvalue -r -df' | " .
		"getinfo -a 'VIEW=$view'";
}
				# no depth buffer, simplest case
exec "@vwraysA | @rtraceA -ffc @res $oct | getinfo -a 'VIEW=$view'";
