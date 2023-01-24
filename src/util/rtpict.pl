#!/usr/bin/perl -w
# RCSid $Id: rtpict.pl,v 2.26 2023/01/24 21:54:49 greg Exp $
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
my %rtraceC = ('-dt',1, '-dc',1, '-dj',1, '-ds',1, '-dr',1, '-dp',1,
		'-ss',1, '-st',1, '-e',1, '-am',1, '-P',1, '-PP',1,
		'-ab',1, '-af',1, '-ai',1, '-aI',1, '-ae',1, '-aE',1,
		'-av',3, '-aw',1, '-aa',1, '-ar',1, '-ad',1, '-as',1,
		'-me',3, '-ma',3, '-mg',1, '-ms',1, '-lr',1, '-lw',1,
		'-ap',2, '-am',1, '-ac',1, '-aC',1);
# boolean rtrace options
my @boolO = ('-w', '-bv', '-dv', '-i', '-u');
# view options and the associated number of arguments
my %vwraysC = ('-vf',1, '-vtv',0, '-vtl',0, '-vth',0, '-vta',0, '-vts',0, '-vtc',0,
		'-x',1, '-y',1, '-vp',3, '-vd',3, '-vu',3, '-vh',1, '-vv',1,
		'-vo',1, '-va',1, '-vs',1, '-vl',1, '-pa',1, '-pj',1, '-pd',1);
# options we need to silently ignore
my %ignoreC = ('-t',1, '-ps',1, '-pt',1, '-pm',1,);
# Starting options for rtrace (rpict values)
my @rtraceA = split(' ', 'rtrace -u- -dt .05 -dc .5 -ds .25 -dr 1 ' .
				'-aa .2 -ar 64 -ad 512 -as 128 -lr 7 -lw 1e-04');
my @vwraysA = ('vwrays', '-pj', '.67');
my @vwrightA = ('vwright', '-vtv');
my @rpictA = ('rpict', '-ps', '1');
my $outpatt = '^-o[vrxlLRXnNsmM]+';
my $refDepth = "";
my $irrad = 0;
my $persist = 0;
my $ambounce = 0;
my $ambcache = 1;
my $ambfile;
my $outlyr;
my $outdir;
my $outpic;
my $outzbf;
OPTION:				# sort through options
while ($#ARGV >= 0 && "$ARGV[0]" =~ /^[-\@]/) {
	# Check for file inclusion
	if ("$ARGV[0]" =~ /^\@/) {
		open my $handle, '<', substr($ARGV[0], 1) or die "No file for $ARGV[0]\n";
		shift @ARGV;
		chomp(my @args = <$handle>);
		close $handle;
		unshift @ARGV, split(/\s+/, "@args");
		next OPTION;
	}
	# Check booleans
	for my $boopt (@boolO) {
		if ("$ARGV[0]" =~ ('^' . $boopt . '[-+01tfynTFYN]?$')) {
			if ("$ARGV[0]" eq '-i') {
				$irrad = ! $irrad;
			} elsif ("$ARGV[0]" =~ /^-i[-+01tfynTFYN]/) {
				$irrad = ("$ARGV[0]" =~ /^-i[+1tyTY]/);
			}
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
		if ("$ARGV[0]" =~ /^-PP?$/) {
			$persist = 1;
		} elsif ("$ARGV[0]" eq '-ab') {
			$ambounce = $ARGV[1];
		} elsif ("$ARGV[0]" eq '-aa') {
			$ambcache = ($ARGV[1] > 0.0);
		} elsif ("$ARGV[0]" eq '-af') {
			$ambfile = "$ARGV[1]";
		}
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
	if ("$ARGV[0]" eq '-o') {
		shift @ARGV;
		$outpic = shift(@ARGV);
	} elsif ("$ARGV[0]" eq '-z') {
		push @rpictA, shift(@ARGV);
		$outzbf = $ARGV[0];
		push @rpictA, shift(@ARGV);
	} elsif ("$ARGV[0]" eq '-n') {
		shift @ARGV;
		$nprocs = shift(@ARGV);
	} elsif ("$ARGV[0]" eq '-d') {
		shift @ARGV;
		$refDepth = ' -d ' . shift(@ARGV);
	} elsif ("$ARGV[0]" =~ "$outpatt") {
		push @rtraceA, $ARGV[0];
		$outlyr = substr(shift(@ARGV), 2);
		$outdir = shift(@ARGV);
	} else {
		die "Unsupported option: " . $ARGV[0] . "\n";
	}
}
die "Number of processes must be positive" if ($nprocs <= 0);
if (defined $outdir) {		# check conflicting options
	die "Options -o and -o* are mutually exclusive\n" if (defined $outpic);
	die "Options -z and -o* are mutually exclusive\n" if (defined $outzbf);
}
if (defined $outpic) {		# redirect output?
	die "File '$outpic' already exists\n" if (-e $outpic);
	open STDOUT, '>', "$outpic";
}
#####################################################################
##### May as well run rpict?
if ($nprocs == 1 && $persist == 0 && !defined($outdir)) {
	push(@rpictA, $ARGV[0]) if ($#ARGV == 0);
	exec @rpictA ;
}
# OK, we're running rtrace in some capacity...
push @rtraceA, ('-n', "$nprocs");
die "Need single octree argument\n" if ($#ARGV != 0);
my $oct = $ARGV[0];
my $view = `@vwrightA 0`;
chomp $view;
my @res = split(/\s/, `@vwraysA -d`);
#####################################################################
##### Resort pixels to reduce ambient cache collisions?
if ($nprocs > 1 && $ambounce > 0 && $ambcache && defined($ambfile)) {
	if (!defined($outzbf) && !defined($outdir)) {
		# Straight picture output, so just shuffle sample order
		system "cnt -s $res[1] $res[3] > /tmp/ord$$.txt";
		die "cnt error\n" if ( $? );
		system "@vwraysA -ff -i < /tmp/ord$$.txt " .
			"| @rtraceA -ffa -ov '$oct' > /tmp/pix$$.txt";
		die "Error running rtrace\n" if ( $? );
		system "( getinfo < /tmp/pix$$.txt | getinfo -a 'VIEW=$view'; " .
			"getinfo - < /tmp/pix$$.txt | rlam /tmp/ord$$.txt - " .
			"| sort -k2rn -k1n ) | pvalue -r -Y $res[3] +X $res[1]";
		die "rlam error\n" if ( $? );
		unlink ("/tmp/ord$$.txt", "/tmp/pix$$.txt");
		exit 0;
	}
	# Else randomize overture calculation to prime ambient cache
	my @ores = (int($res[1]/6), int($res[3]/6));
	print STDERR "Running $ores[0] by $ores[1] overture calculation " .
			"to populate '$ambfile'...\n";
	system "cnt -s @ores | @vwraysA -i -ff -x $ores[0] -y $ores[1] -pj 0 " .
		"| @rtraceA -ff -ov '$oct' > /dev/null";
	die "Failure running overture\n" if ( $? );
	print STDERR "Finished overture.\n";
}
#####################################################################
##### Generating picture with depth buffer?
if (defined $outzbf) {
	exec "@vwraysA -ff | @rtraceA -fff -olv @res '$oct' | " .
		"rsplit -ih -iH -f -of '$outzbf' -oh -oH -of3 - | " .
		"pvalue -r -df | getinfo -a 'VIEW=$view'";
}
#####################################################################
##### Base case with output picture only?
if (! defined $outdir) {
	exec "@vwraysA -ff | @rtraceA -ffc @res '$oct' | getinfo -a 'VIEW=$view'";
}
#####################################################################
##### Layered image output case
my @rsplitA = ("rsplit", "'-t	'", "-ih", "-iH", "-oh", "-oH");
# Supported rtrace -o* options and output file name.typ
my %rtoutC = (
	v =>	'radiance.hdr',
	r =>	'r_refl.hdr',
	x =>	'r_unrefl.hdr',
	l =>	'd_effective.dpt',
	L =>	'd_firstsurf.dpt',
	R =>	'd_refl.dpt',
	X =>	'd_unrefl.dpt',
	n =>	'perturbed.nrm',
	N =>	'unperturbed.nrm',
	s =>	'surface.idx',
	m =>	'modifier.idx',
	M =>	'material.idx'
);
# Arguments for rsplit based on output file type
my %rcodeC = (
	'.hdr',	['-of3', '!pvalue -r -df -u'],
	'.dpt',	['-of', "!rcode_depth$refDepth -ff"],
	'.nrm',	['-of3', '!rcode_norm -ff'],
	'.idx',	['-oa', '!rcode_ident "-t	"']
);
if ($irrad) {		# adjust actions for -i+ option
	$rtoutC{'v'} = 'irradiance.hdr';
	delete $rtoutC{'r'};
	delete $rtoutC{'x'};
	delete $rtoutC{'R'};
	delete $rtoutC{'X'};
}
if (! -d $outdir) {
	mkdir $outdir or die("Cannot create output directory '$outdir'\n");
}
			# construct rsplit command
foreach my $oval (split //, $outlyr) {
	die "Duplicate or unsupported type '$oval' in -o$outlyr\n"
					if (!defined $rtoutC{$oval});
	my $outfile = "$outdir/$rtoutC{$oval}";
	die "File '$outfile' already exists\n" if (-e $outfile);
	my ($otyp) = ($rtoutC{$oval} =~ /(\.[^.]+)$/);
	push @rsplitA, $rcodeC{$otyp}[0];
	push @rsplitA, qq{'$rcodeC{$otyp}[1] > "$outfile"'};
	delete $rtoutC{$oval};
}
			# call rtrace + rsplit
exec "@vwraysA -ff | @rtraceA -fff @res '$oct' | getinfo -a 'VIEW=$view' | @rsplitA";
