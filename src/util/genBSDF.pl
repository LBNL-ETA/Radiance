#!/usr/bin/perl -w
# RCSid $Id: genBSDF.pl,v 2.18 2011/06/01 22:50:24 greg Exp $
#
# Compute BSDF based on geometry and material description
#
#	G. Ward
#
use strict;
use File::Temp qw/ :mktemp  /;
sub userror {
	print STDERR "Usage: genBSDF [-n Nproc][-c Nsamp][-t{3|4} Nlog2][-r \"ropts\"][-dim xmin xmax ymin ymax zmin zmax][{+|-}f][{+|-}b][{+|-}mgf][{+|-}geom] [input ..]\n";
	exit 1;
}
my $td = mkdtemp("/tmp/genBSDF.XXXXXX");
chomp $td;
my $tensortree = 0;
my $ttlog2 = 4;
my $nsamp = 1000;
my $rtargs = "-w -ab 5 -ad 700 -lw 3e-6";
my $mgfin = 0;
my $geout = 1;
my $nproc = 1;
my $doforw = 0;
my $doback = 1;
my @dim;
# Get options
while ($#ARGV >= 0) {
	if ("$ARGV[0]" =~ /^[-+]m/) {
		$mgfin = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" eq "-r") {
		$rtargs = "$rtargs $ARGV[1]";
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^[-+]g/) {
		$geout = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]f/) {
		$doforw = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]b/) {
		$doback = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^-t[34]$/) {
		$tensortree = substr($ARGV[0], 2, 1);
		$ttlog2 = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-c") {
		$nsamp = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-n") {
		$nproc = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-d/) {
		userror() if ($#ARGV < 6);
		@dim = @ARGV[1..6];
		shift @ARGV for (1..6);
	} elsif ("$ARGV[0]" =~ /^[-+]./) {
		userror();
	} else {
		last;
	}
	shift @ARGV;
}
# Check that we're actually being asked to do something
die "Must have at least one of +forward or +backward\n" if (!$doforw && !$doback);
# Name our own persist file?
my $persistfile;
if ( $tensortree && $nproc > 1 && "$rtargs" !~ /-PP /) {
	$persistfile = "$td/pfile.txt";
	$rtargs = "-PP $persistfile $rtargs";
}
# Get scene description and dimensions
my $radscn = "$td/device.rad";
my $mgfscn = "$td/device.mgf";
my $octree = "$td/device.oct";
if ( $mgfin ) {
	system "mgfilt '#,o,xf,c,cxy,cspec,cmix,m,sides,rd,td,rs,ts,ir,v,p,n,f,fh,sph,cyl,cone,prism,ring,torus' @ARGV > $mgfscn";
	die "Could not load MGF input\n" if ( $? );
	system "mgf2rad $mgfscn > $radscn";
} else {
	system "xform -e @ARGV > $radscn";
	die "Could not load Radiance input\n" if ( $? );
	system "rad2mgf $radscn > $mgfscn" if ( $geout );
}
if ($#dim != 5) {
	@dim = split ' ', `getbbox -h $radscn`;
}
print STDERR "Warning: Device extends into room\n" if ($dim[5] > 1e-5);
# Add receiver surfaces (rectangular)
my $fmodnm="receiver_face";
my $bmodnm="receiver_behind";
open(RADSCN, ">> $radscn");
print RADSCN "void glow $fmodnm\n0\n0\n4 1 1 1 0\n\n";
print RADSCN "$fmodnm source f_receiver\n0\n0\n4 0 0 1 180\n";
print RADSCN "void glow $bmodnm\n0\n0\n4 1 1 1 0\n\n";
print RADSCN "$bmodnm source b_receiver\n0\n0\n4 0 0 -1 180\n";
close RADSCN;
# Generate octree
system "oconv -w $radscn > $octree";
die "Could not compile scene\n" if ( $? );
# Output XML prologue
print
'<?xml version="1.0" encoding="UTF-8"?>
<WindowElement xmlns="http://windows.lbl.gov" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://windows.lbl.gov/BSDF-v1.4.xsd">
<WindowElementType>System</WindowElementType>
<Optical>
<Layer>
	<Material>
		<Name>Name</Name>
		<Manufacturer>Manufacturer</Manufacturer>
';
printf "\t\t<Thickness unit=\"Meter\">%.3f</Thickness>\n", $dim[5] - $dim[4];
printf "\t\t<Width unit=\"Meter\">%.3f</Width>\n", $dim[1] - $dim[0];
printf "\t\t<Height unit=\"Meter\">%.3f</Height>\n", $dim[3] - $dim[2];
print "\t\t<DeviceType>Integral</DeviceType>\n";
# Output MGF description if requested
if ( $geout ) {
	print "\t\t<Geometry format=\"MGF\" unit=\"Meter\">\n";
	printf "xf -t %.6f %.6f 0\n", -($dim[0]+$dim[1])/2, -($dim[2]+$dim[3])/2;
	open(MGFSCN, "< $mgfscn");
	while (<MGFSCN>) { print $_; }
	close MGFSCN;
	print "xf\n";
	print "\t\t</Geometry>\n";
}
print "	</Material>\n";
# Set up surface sampling
my $nx = int(sqrt($nsamp*($dim[1]-$dim[0])/($dim[3]-$dim[2])) + .5);
my $ny = int($nsamp/$nx + .5);
$nsamp = $nx * $ny;
my $ns = 2**$ttlog2;
my (@pdiv, $disk2sq, $sq2disk, $tcal, $kcal);
# Create data segments (all the work happens here)
if ( $tensortree ) {
	do_tree_bsdf();
} else {
	do_matrix_bsdf();
}
# Output XML epilogue
print
'</Layer>
</Optical>
</WindowElement>
';
# Clean up temporary files and exit
if ( $persistfile && open(PFI, "< $persistfile") ) {
	while (<PFI>) {
		s/^[^ ]* //;
		kill('ALRM', $_);
		last;
	}
	close PFI;
}
exec("rm -rf $td");

#-------------- End of main program segment --------------#

#++++++++++++++ Tensor tree BSDF generation ++++++++++++++#
sub do_tree_bsdf {
# Get sampling rate and subdivide task
my $ns2 = $ns;
$ns2 /= 2 if ( $tensortree == 3 );
@pdiv = (0, int($ns2/$nproc));
my $nrem = $ns2 % $nproc;
for (my $i = 1; $i < $nproc; $i++) {
	my $nv = $pdiv[$i] + $pdiv[1];
	++$nv if ( $nrem-- > 0 );
	push @pdiv, $nv;
}
die "Script error 1" if ($pdiv[-1] != $ns2);
# Shirley-Chiu mapping from unit square to disk
$sq2disk = '
in_square_a = 2*in_square_x - 1;
in_square_b = 2*in_square_y - 1;
in_square_rgn = if(in_square_a + in_square_b,
			if(in_square_a - in_square_b, 1, 2),
			if(in_square_b - in_square_a, 3, 4));
out_disk_r = .999995*select(in_square_rgn, in_square_a, in_square_b,
			-in_square_a, -in_square_b);
out_disk_phi = PI/4 * select(in_square_rgn,
				in_square_b/in_square_a,
				2 - in_square_a/in_square_b,
				4 + in_square_b/in_square_a,
				if(in_square_b*in_square_b,
					6 - in_square_a/in_square_b, 0));
Dx = out_disk_r*cos(out_disk_phi);
Dy = out_disk_r*sin(out_disk_phi);
Dz = sqrt(1 - out_disk_r*out_disk_r);
';
# Shirley-Chiu mapping from unit disk to square
$disk2sq = '
norm_radians(p) : if(-p - PI/4, p + 2*PI, p);
in_disk_r = .999995*sqrt(Dx*Dx + Dy*Dy);
in_disk_phi = norm_radians(atan2(Dy, Dx));
in_disk_rgn = floor((in_disk_phi + PI/4)/(PI/2)) + 1;
out_square_a = select(in_disk_rgn,
			in_disk_r,
			(PI/2 - in_disk_phi)*in_disk_r/(PI/4),
			-in_disk_r,
			(in_disk_phi - 3*PI/2)*in_disk_r/(PI/4));
out_square_b = select(in_disk_rgn,
			in_disk_phi*in_disk_r/(PI/4),
			in_disk_r,
			(PI - in_disk_phi)*in_disk_r/(PI/4),
			-in_disk_r);
out_square_x = (out_square_a + 1)/2;
out_square_y = (out_square_b + 1)/2;
';
# Announce ourselves in XML output
print "\t<DataDefinition>\n";
print "\t\t<IncidentDataStructure>TensorTree$tensortree</IncidentDataStructure>\n";
print "\t</DataDefinition>\n";
# Fork parallel rtcontrib processes to compute each side
if ( $doback ) {
	for (my $proc = 0; $proc < $nproc; $proc++) {
		bg_tree_rtcontrib(0, $proc);
	}
	while (wait() >= 0) {
		die "rtcontrib process reported error" if ( $? );
	}
	ttree_out(0);
}
if ( $doforw ) {
	for (my $proc = 0; $proc < $nproc; $proc++) {
		bg_tree_rtcontrib(1, $proc);
	}
	while (wait() >= 0) {
		die "rtcontrib process reported error" if ( $? );
	}
	ttree_out(1);
}
}	# end of sub do_tree_bsdf()

# Run i'th rtcontrib process for generating tensor tree samples
sub bg_tree_rtcontrib {
	my $pid = fork();
	die "Cannot fork new process" unless defined $pid;
	if ($pid > 0) { return $pid; }
	my $forw = shift;
	my $pn = shift;
	my $pbeg = $pdiv[$pn];
	my $plen = $pdiv[$pn+1] - $pbeg;
	my $matargs = "-m $bmodnm";
	if ( !$forw || !$doback ) { $matargs .= " -m $fmodnm"; }
	my $cmd = "rtcontrib $rtargs -h -ff -fo -c $nsamp " .
		"-e '$disk2sq' -bn '$ns*$ns' " .
		"-b '$ns*floor(out_square_x*$ns)+floor(out_square_y*$ns)' " .
		"-o $td/%s_" . sprintf("%03d", $pn) . ".flt $matargs $octree";
	if ( $tensortree == 3 ) {
		# Isotropic BSDF
		$cmd = "cnt $plen $ny $nx " .
			"| rcalc -e 'r1=rand(($pn+.8681)*recno-.673892)' " .
			"-e 'r2=rand(($pn-5.37138)*recno+67.1737811)' " .
			"-e 'r3=rand(($pn+3.17603772)*recno+83.766771)' " .
			"-e 'Dx=1-2*($pbeg+\$1+r1)/$ns;Dy:0;Dz=sqrt(1-Dx*Dx)' " .
			"-e 'xp=(\$3+r2)*(($dim[1]-$dim[0])/$nx)+$dim[0]' " .
			"-e 'yp=(\$2+r3)*(($dim[3]-$dim[2])/$ny)+$dim[2]' " .
			"-e 'zp=$dim[5-$forw]' -e 'myDz=Dz*($forw*2-1)' " .
			"-e '\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz' " .
			"-e '\$4=Dx;\$5=Dy;\$6=myDz' -of " .
			"| $cmd";
	} else {
		# Anisotropic BSDF
		# Sample area vertically to improve load balance, since
		# shading systems usually have bilateral symmetry (L-R)
		$cmd = "cnt $plen $ns $ny $nx " .
			"| rcalc -e 'r1=rand(($pn+.8681)*recno-.673892)' " .
			"-e 'r2=rand(($pn-5.37138)*recno+67.1737811)' " .
			"-e 'r3=rand(($pn+3.17603772)*recno+83.766771)' " .
			"-e 'r4=rand(($pn-2.3857833)*recno-964.72738)' " .
			"-e 'in_square_x=($pbeg+\$1+r1)/$ns' " .
			"-e 'in_square_y=(\$2+r2)/$ns' -e '$sq2disk' " .
			"-e 'xp=(\$4+r3)*(($dim[1]-$dim[0])/$nx)+$dim[0]' " .
			"-e 'yp=(\$3+r4)*(($dim[3]-$dim[2])/$ny)+$dim[2]' " .
			"-e 'zp=$dim[5-$forw]' -e 'myDz=Dz*($forw*2-1)' " .
			"-e '\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz' " .
			"-e '\$4=Dx;\$5=Dy;\$6=myDz' -of " .
			"| $cmd";
	}
# print STDERR "Starting: $cmd\n";
	exec($cmd);		# no return; status report to parent via wait
	die "Cannot exec: $cmd\n";
}	# end of bg_tree_rtcontrib()

# Simplify and output tensor tree results
sub ttree_out {
	my $forw = shift;
	my $side = ("Back","Front")[$forw];
# Only output one transmitted distribution, preferring backwards
if ( !$forw || !$doback ) {
print
'	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
			<WavelengthDataDirection>Transmission</WavelengthDataDirection>
			<AngleBasis>LBNL/Shirley-Chiu</AngleBasis>
			<ScatteringDataType>BTDF</ScatteringDataType>
			<ScatteringData>
';
system "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
	q{-e '$1=(0.265*$1+0.670*$2+0.065*$3)/Omega' -of } .
	"$td/" . ($bmodnm,$fmodnm)[$forw] . "_???.flt " .
	"| rttree_reduce -h -ff -r $tensortree -g $ttlog2";
die "Failure running rttree_reduce" if ( $? );
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
';
}
# Output reflection
print
'	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
';
print "\t\t\t<WavelengthDataDirection>Reflection $side</WavelengthDataDirection>\n";
print '
			<AngleBasis>LBNL/Shirley-Chiu</AngleBasis>
			<ScatteringDataType>BRDF</ScatteringDataType>
			<ScatteringData>
';
system "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
	q{-e '$1=(0.265*$1+0.670*$2+0.065*$3)/Omega' -of } .
	"$td/" . ($fmodnm,$bmodnm)[$forw] . "_???.flt " .
	"| rttree_reduce -h -ff -r $tensortree -g $ttlog2";
die "Failure running rttree_reduce" if ( $? );
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
';
}	# end of ttree_out()

#------------- End of do_tree_bsdf() & subroutines -------------#

#+++++++++++++++ Klems matrix BSDF generation +++++++++++++++#
sub do_matrix_bsdf {
# Set up sampling of portal
# Kbin to produce incident direction in full Klems basis with (x1,x2) randoms
$tcal = '
DEGREE : PI/180;
sq(x) : x*x;
Kpola(r) : select(r+1, -5, 5, 15, 25, 35, 45, 55, 65, 75, 90);
Knaz(r) : select(r, 1, 8, 16, 20, 24, 24, 24, 16, 12);
Kaccum(r) : if(r-.5, Knaz(r) + Kaccum(r-1), 0);
Kmax : Kaccum(Knaz(0));
Kfindrow(r, rem) : if(rem-Knaz(r)+.5, Kfindrow(r+1, rem-Knaz(r)), r);
Krow = if(Kbin-(Kmax-.5), 0, Kfindrow(1, Kbin));
Kcol = Kbin - Kaccum(Krow-1);
Kazi = 360*DEGREE * (Kcol + (.5 - x2)) / Knaz(Krow);
Kpol = DEGREE * (x1*Kpola(Krow) + (1-x1)*Kpola(Krow-1));
sin_kpol = sin(Kpol);
Dx = cos(Kazi)*sin_kpol;
Dy = sin(Kazi)*sin_kpol;
Dz = sqrt(1 - sin_kpol*sin_kpol);
KprojOmega = PI * if(Kbin-.5,
	(sq(cos(Kpola(Krow-1)*DEGREE)) - sq(cos(Kpola(Krow)*DEGREE)))/Knaz(Krow),
	1 - sq(cos(Kpola(1)*DEGREE)));
';
# Compute Klems bin from exiting ray direction (forward or backward)
$kcal = '
DEGREE : PI/180;
abs(x) : if(x, x, -x);
Acos(x) : 1/DEGREE * if(x-1, 0, if(-1-x, 0, acos(x)));
posangle(a) : if(-a, a + 2*PI, a);
Atan2(y,x) : 1/DEGREE * posangle(atan2(y,x));
kpola(r) : select(r, 5, 15, 25, 35, 45, 55, 65, 75, 90);
knaz(r) : select(r, 1, 8, 16, 20, 24, 24, 24, 16, 12);
kaccum(r) : if(r-.5, knaz(r) + kaccum(r-1), 0);
kfindrow(r, pol) : if(r-kpola(0)+.5, r,
		if(pol-kpola(r), kfindrow(r+1, pol), r) );
kazn(azi,inc) : if((360-.5*inc)-azi, floor((azi+.5*inc)/inc), 0);
kbin2(pol,azi) = select(kfindrow(1, pol),
		kazn(azi,360/knaz(1)),
		kaccum(1) + kazn(azi,360/knaz(2)),
		kaccum(2) + kazn(azi,360/knaz(3)),
		kaccum(3) + kazn(azi,360/knaz(4)),
		kaccum(4) + kazn(azi,360/knaz(5)),
		kaccum(5) + kazn(azi,360/knaz(6)),
		kaccum(6) + kazn(azi,360/knaz(7)),
		kaccum(7) + kazn(azi,360/knaz(8)),
		kaccum(8) + kazn(azi,360/knaz(9))
	);
kbin = kbin2(Acos(abs(Dz)),Atan2(Dy,Dx));
';
my $ndiv = 145;
# Compute scattering data using rtcontrib
my @tfarr;
my @rfarr;
my @tbarr;
my @rbarr;
my $cmd;
my $rtcmd = "rtcontrib $rtargs -h -ff -fo -n $nproc -c $nsamp " .
	"-e '$kcal' -b kbin -bn $ndiv " .
	"-o '$td/%s.flt' -m $fmodnm -m $bmodnm $octree";
my $rccmd = "rcalc -e '$tcal' " .
	"-e 'mod(n,d):n-floor(n/d)*d' -e 'Kbin=mod(recno-.999,$ndiv)' " .
	q{-if3 -e '$1=(0.265*$1+0.670*$2+0.065*$3)/KprojOmega'};
if ( $doforw ) {
$cmd = "cnt $ndiv $ny $nx | rcalc -of -e '$tcal' " .
	"-e 'xp=(\$3+rand(.12*recno+288))*(($dim[1]-$dim[0])/$nx)+$dim[0]' " .
	"-e 'yp=(\$2+rand(.37*recno-44))*(($dim[3]-$dim[2])/$ny)+$dim[2]' " .
	"-e 'zp:$dim[4]' " .
	q{-e 'Kbin=$1;x1=rand(2.75*recno+3.1);x2=rand(-2.01*recno-3.37)' } .
	q{-e '$1=xp-Dx;$2=yp-Dy;$3=zp-Dz;$4=Dx;$5=Dy;$6=Dz' } .
	"| $rtcmd";
system "$cmd" || die "Failure running: $cmd\n";
@tfarr = `$rccmd $td/$fmodnm.flt`;
die "Failure running: $rccmd $td/$fmodnm.flt\n" if ( $? );
@rfarr = `$rccmd $td/$bmodnm.flt`;
die "Failure running: $rccmd $td/$bmodnm.flt\n" if ( $? );
}
if ( $doback ) {
$cmd = "cnt $ndiv $ny $nx | rcalc -of -e '$tcal' " .
	"-e 'xp=(\$3+rand(.35*recno-15))*(($dim[1]-$dim[0])/$nx)+$dim[0]' " .
	"-e 'yp=(\$2+rand(.86*recno+11))*(($dim[3]-$dim[2])/$ny)+$dim[2]' " .
	"-e 'zp:$dim[5]' " .
	q{-e 'Kbin=$1;x1=rand(1.21*recno+2.75);x2=rand(-3.55*recno-7.57)' } .
	q{-e '$1=xp-Dx;$2=yp-Dy;$3=zp+Dz;$4=Dx;$5=Dy;$6=-Dz' } .
	"| $rtcmd";
system "$cmd" || die "Failure running: $cmd\n";
@tbarr = `$rccmd $td/$bmodnm.flt`;
die "Failure running: $rccmd $td/$bmodnm.flt\n" if ( $? );
@rbarr = `$rccmd $td/$fmodnm.flt`;
die "Failure running: $rccmd $td/$fmodnm.flt\n" if ( $? );
}
# Output angle basis
print
'	<DataDefinition>
		<IncidentDataStructure>Columns</IncidentDataStructure>
		<AngleBasis>
			<AngleBasisName>LBNL/Klems Full</AngleBasisName>
			<AngleBasisBlock>
				<Theta>0</Theta>
				<nPhis>1</nPhis>
				<ThetaBounds>
					<LowerTheta>0</LowerTheta>
					<UpperTheta>5</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>10</Theta>
				<nPhis>8</nPhis>
				<ThetaBounds>
					<LowerTheta>5</LowerTheta>
					<UpperTheta>15</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>20</Theta>
				<nPhis>16</nPhis>
				<ThetaBounds>
					<LowerTheta>15</LowerTheta>
					<UpperTheta>25</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>30</Theta>
				<nPhis>20</nPhis>
				<ThetaBounds>
					<LowerTheta>25</LowerTheta>
					<UpperTheta>35</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>40</Theta>
				<nPhis>24</nPhis>
				<ThetaBounds>
					<LowerTheta>35</LowerTheta>
					<UpperTheta>45</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>50</Theta>
				<nPhis>24</nPhis>
				<ThetaBounds>
					<LowerTheta>45</LowerTheta>
					<UpperTheta>55</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>60</Theta>
				<nPhis>24</nPhis>
				<ThetaBounds>
					<LowerTheta>55</LowerTheta>
					<UpperTheta>65</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>70</Theta>
				<nPhis>16</nPhis>
				<ThetaBounds>
					<LowerTheta>65</LowerTheta>
					<UpperTheta>75</UpperTheta>
				</ThetaBounds>
				</AngleBasisBlock>
				<AngleBasisBlock>
				<Theta>82.5</Theta>
				<nPhis>12</nPhis>
				<ThetaBounds>
					<LowerTheta>75</LowerTheta>
					<UpperTheta>90</UpperTheta>
				</ThetaBounds>
			</AngleBasisBlock>
		</AngleBasis>
	</DataDefinition>
';
if ( $doforw ) {
print
'	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
			<WavelengthDataDirection>Transmission Front</WavelengthDataDirection>
			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>
			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>
			<ScatteringDataType>BTDF</ScatteringDataType>
			<ScatteringData>
';
# Output front transmission (transposed order)
for (my $od = 0; $od < $ndiv; $od++) {
	for (my $id = 0; $id < $ndiv; $id++) {
		print $tfarr[$ndiv*$id + $od];
	}
	print "\n";
}
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
			<WavelengthDataDirection>Reflection Front</WavelengthDataDirection>
			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>
			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>
			<ScatteringDataType>BRDF</ScatteringDataType>
			<ScatteringData>
';
# Output front reflection (transposed order)
for (my $od = 0; $od < $ndiv; $od++) {
	for (my $id = 0; $id < $ndiv; $id++) {
		print $rfarr[$ndiv*$id + $od];
	}
	print "\n";
}
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
';
}
if ( $doback ) {
print
'	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
			<WavelengthDataDirection>Transmission Back</WavelengthDataDirection>
			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>
			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>
			<ScatteringDataType>BTDF</ScatteringDataType>
			<ScatteringData>
';
# Output back transmission (transposed order)
for (my $od = 0; $od < $ndiv; $od++) {
	for (my $id = 0; $id < $ndiv; $id++) {
		print $tbarr[$ndiv*$id + $od];
	}
	print "\n";
}
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
	<WavelengthData>
		<LayerNumber>System</LayerNumber>
		<Wavelength unit="Integral">Visible</Wavelength>
		<SourceSpectrum>CIE Illuminant D65 1nm.ssp</SourceSpectrum>
		<DetectorSpectrum>ASTM E308 1931 Y.dsp</DetectorSpectrum>
		<WavelengthDataBlock>
			<WavelengthDataDirection>Reflection Back</WavelengthDataDirection>
			<ColumnAngleBasis>LBNL/Klems Full</ColumnAngleBasis>
			<RowAngleBasis>LBNL/Klems Full</RowAngleBasis>
			<ScatteringDataType>BRDF</ScatteringDataType>
			<ScatteringData>
';
# Output back reflection (transposed order)
for (my $od = 0; $od < $ndiv; $od++) {
	for (my $id = 0; $id < $ndiv; $id++) {
		print $rbarr[$ndiv*$id + $od];
	}
	print "\n";
}
print
'			</ScatteringData>
		</WavelengthDataBlock>
	</WavelengthData>
';
}
}
#------------- End of do_matrix_bsdf() --------------#
