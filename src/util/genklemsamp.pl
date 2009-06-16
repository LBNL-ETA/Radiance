#!/usr/bin/perl -w
# RCSid $Id: genklemsamp.pl,v 2.2 2009/06/16 04:54:08 greg Exp $
#
# Sample Klems (full) directions impinging on surface(s)
#
#	G. Ward
#
if ($#ARGV < 0) {
	print STDERR, "Usage: genklemsamp [-c N ][-f{a|f|d}] [view opts] [geom.rad ..]\n";
	exit 1;
}
my $td = `mktemp -d /tmp/genklemsamp.XXXXXX`;
chomp $td;
my $nsamp = 1000;
my $fmt = "a";
my @vopts="-vo 1e-5";
while ($#ARGV >= 0) {
	if ("$ARGV[0]" eq "-c") {
		shift @ARGV;
		$nsamp = $ARGV[0];
	} elsif ("$ARGV[0]" =~ /^-f[afd]$/) {
		$fmt = substr($ARGV[0], 2, 1);
	} elsif ("$ARGV[0]" =~ /^-v[pdu]$/) {
		push @vopts, "@ARGV[0..3]";
		shift @ARGV; shift @ARGV; shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-v[fhvo]$/) {
		push @vopts, "@ARGV[0..1]";
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-v./) {
		print STDERR, "Unsupported view option: $ARGV[0]\n";
		exit 1;
	} elsif ("$ARGV[0]" =~ /^-./) {
		print STDERR, "Unknown option: $ARGV[0]\n";
		exit 1;
	} else {
		last;
	}
	shift @ARGV;
}
my $ofmt = "";
$ofmt = "-o$fmt" if ("$fmt" ne "a");
# Require parallel view and compile settings
push @vopts, "-vtl";
my $vwset = `vwright @vopts V`;
$_ = $vwset; s/^.*Vhn://; s/;.*$//;
my $width = $_;
$_ = $vwset; s/^.*Vvn://; s/;.*$//;
my $height = $_;
my $sca = sqrt($nsamp/($width*$height));
# Kbin is input to get direction in full Klems basis with (x1,x2) randoms
my $tcal = '
DEGREE : PI/180;
Kpola(r) : select(r+1, -5, 5, 15, 25, 35, 45, 55, 65, 75, 90);
Knaz(r) : select(r, 1, 8, 16, 20, 24, 24, 24, 16, 12);
Kaccum(r) : if(r-.5, Knaz(r) + Kaccum(r-1), 0);
Kmax : Kaccum(Knaz(0));
Kfindrow(r, rem) : if(rem-Knaz(r)-.5, Kfindrow(r+1, rem-Knaz(r)), r);
Krow = if(Kbin-(Kmax+.5), 0, Kfindrow(1, Kbin));
Kcol = Kbin - Kaccum(Krow-1) - 1;
Kazi = 360*DEGREE * (Kcol + .5 - x2) / Knaz(Krow);
Kpol = DEGREE * (x1*Kpola(Krow) + (1-x1)*Kpola(Krow-1));
sin_kpol = sin(Kpol);
K0 = cos(Kazi)*sin_kpol;
K1 = sin(Kazi)*sin_kpol;
K2 = cos(Kpol);
';
my $ndiv = 145;
# Do we have any Radiance input files?
if ($#ARGV >= 0) {
	system "oconv -f @ARGV > $td/surf.oct";
	# Set our own view center and size based on bounding cube
	$_ = $vwset; s/^.*Vdx://; s/;.*$//;
	my @vd = $_;
	$_ = $vwset; s/^.*Vdy://; s/;.*$//;
	push @vd, $_;
	$_ = $vwset; s/^.*Vdz://; s/;.*$//;
	push @vd, $_;
	my @bcube = split /\s+/, `getinfo -d < $td/surf.oct`;
	$width = $bcube[3]*sqrt(3);
	$height = $width;
	push @vopts, ("-vp", $bcube[0]+$bcube[3]/2-$width/2*$vd[0],
			$bcube[1]+$bcube[3]/2-$width/2*$vd[1],
			$bcube[2]+$bcube[3]/2-$width/2*$vd[2]);
	push @vopts, ("-vh", $width, "-vh", $height);
	$vwset = `vwright @vopts V`;
	my $xres;
	my $yres;
	my $ntot;
	# This generally passes through the loop twice to get density right
	do {
		$xres = int($width*$sca) + 1;
		$yres = int($height*$sca) + 1;
		system "vwrays -ff -x $xres -y $yres -pa 0 -pj .7 @vopts " .
			"| rtrace -w -h -ff -opLN $td/surf.oct " .
			q{| rcalc -e 'cond=and(5e9-$4,nOK);$1=$1;$2=$2;$3=$3' } .
				"-e 'and(a,b):if(a,b,a);sq(x):x*x' -e '$vwset' " .
				q{-e 'nOK=sq(Vdx*$5+Vdy*$6+Vdz*$7)-.999' } .
				"-if7 -of > $td/origins.flt";
		$ntot = -s "$td/origins.flt";
		$ntot /= 3*4;
		if ($ntot == 0) {
			print STDERR, "View direction does not correspond to any surfaces\n";
			exit 1;
		}
		$sca *= 1.05 * sqrt($nsamp/$ntot);
	} while ($ntot < $nsamp);
	# All set to produce our samples
	for ($k = 1; $k <= $ndiv; $k++) {
		my $rn = rand(10);
		my $r1 = rand; my $r2 = rand;
		# Chance of using = (number_still_needed)/(number_left_avail)
		system "rcalc $ofmt -e '$tcal' -e '$vwset' " .
			"-e 'cond=($nsamp-outno+1)/($ntot-recno+1)-rand($rn+recno)' " .
			"-e 'Kbin=$k;x1=rand($r1+recno);x2=rand($r2+recno)' " .
			q{-e '$1=$1+Vo*Vdx; $2=$2+Vo*Vdy; $3=$3+Vo*Vdz' } .
			q{-e '$4=K0*Vhx+K1*Vvx+K2*Vdx' } .
			q{-e '$5=K0*Vhy+K1*Vvy+K2*Vdy' } .
			q{-e '$6=K0*Vhz+K1*Vvz+K2*Vdz' } .
			"-if3 $td/origins.flt";
	}
} else {
	# No Radiance input files, so just sample over parallel view
	my $xres = int($width*$sca) + 1;
	my $yres = int($height*$sca) + 1;
	my $ntot = $xres * $yres;
	for ($k = 1; $k <= $ndiv; $k++) {
		my $rn = rand(10);
		my $r1 = rand; my $r2 = rand;
		my $r3 = rand; my $r4 = rand;
		# Chance of using = (number_still_needed)/(number_left_avail)
		system "cnt $yres $xres | rcalc $ofmt -e '$tcal' -e '$vwset' " .
			q{-e 'xc=$2;yc=$1' } .
			"-e 'cond=($nsamp-outno+1)/($ntot-recno+1)-rand($rn+recno)' " .
			"-e 'Kbin=$k;x1=rand($r1+recno);x2=rand($r2+recno)' " .
			"-e 'hpos=Vhn*((xc+rand($r3+recno))/$xres - .5)' " .
			"-e 'vpos=Vvn*((yc+rand($r4+recno))/$yres - .5)' " .
			q{-e '$1=Vpx+Vo*Vdx+hpos*Vhx+vpos*Vvx' } .
			q{-e '$2=Vpy+Vo*Vdy+hpos*Vhy+vpos*Vvy' } .
			q{-e '$3=Vpz+Vo*Vdz+hpos*Vhz+vpos*Vvz' } .
			q{-e '$4=K0*Vhx+K1*Vvx+K2*Vdx' } .
			q{-e '$5=K0*Vhy+K1*Vvy+K2*Vdy' } .
			q{-e '$6=K0*Vhz+K1*Vvz+K2*Vdz' } ;
	}
}
system "rm -rf $td";
