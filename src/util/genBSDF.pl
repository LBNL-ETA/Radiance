#!/usr/bin/perl -w
# RCSid $Id: genBSDF.pl,v 2.5 2010/12/09 17:00:43 greg Exp $
#
# Compute BSDF based on geometry and material description
#
#	G. Ward
#
use strict;
sub userror {
	print STDERR "Usage: genBSDF [-n Nproc][-c Nsamp][-dim xmin xmax ymin ymax zmin zmax][{+|-}mgf][{+|-}geom] [input ..]\n";
	exit 1;
}
my $td = `mktemp -d /tmp/genBSDF.XXXXXX`;
chomp $td;
my $nsamp = 1000;
my $mgfin = 0;
my $geout = 1;
my $nproc = 1;
my @dim;
# Get options
while ($#ARGV >= 0) {
	if ("$ARGV[0]" =~ /^[-+]m/) {
		$mgfin = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]g/) {
		$geout = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" eq "-c") {
		$nsamp = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-n") {
		$nproc = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-d/) {
		userror() if ($#ARGV < 6);
		@dim = "@ARGV[1..6]";
		shift @ARGV for (1..6);
	} elsif ("$ARGV[0]" =~ /^[-+]./) {
		userror();
	} else {
		last;
	}
	shift @ARGV;
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
	system "cat @ARGV | xform -e > $radscn";
	die "Could not load Radiance input\n" if ( $? );
	system "rad2mgf $radscn > $mgfscn" if ( $geout );
}
if ($#dim != 5) {
	@dim = split /\s+/, `getbbox -h $radscn`;
	shift @dim;
}
print STDERR "Warning: Device extends into room\n" if ($dim[5] > 1e-5);
# Add receiver surface (rectangle)
my $modnm="_receiver_black_";
open(RADSCN, ">> $radscn");
print RADSCN "void glow $modnm\n0\n0\n4 0 0 0 0\n\n";
print RADSCN "$modnm polygon _receiver_\n0\n0\n12\n";
print RADSCN "\t",$dim[0],"\t",$dim[2],"\t",$dim[5]+1e-5,"\n";
print RADSCN "\t",$dim[0],"\t",$dim[3],"\t",$dim[5]+1e-5,"\n";
print RADSCN "\t",$dim[1],"\t",$dim[3],"\t",$dim[5]+1e-5,"\n";
print RADSCN "\t",$dim[1],"\t",$dim[2],"\t",$dim[5]+1e-5,"\n";
close RADSCN;
# Generate octree
system "oconv -w $radscn > $octree";
die "Could not compile scene\n" if ( $? );
# Set up sampling
# Kbin to produce incident direction in full Klems basis with (x1,x2) randoms
my $tcal = '
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
Dx = -cos(Kazi)*sin_kpol;
Dy = sin(Kazi)*sin_kpol;
Dz = sqrt(1 - sin_kpol*sin_kpol);
KprojOmega = PI * if(Kbin-.5,
	(sq(cos(Kpola(Krow-1)*DEGREE)) - sq(cos(Kpola(Krow)*DEGREE)))/Knaz(Krow),
	1 - sq(cos(Kpola(1)*DEGREE)));
';
# Compute Klems bin from exiting ray direction
my $kcal = '
DEGREE : PI/180;
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
kbin = kbin2(Acos(Dz), Atan2(Dy, -Dx));
';
my $ndiv = 145;
my $nx = int(sqrt($nsamp*($dim[1]-$dim[0])/($dim[3]-$dim[2])) + .5);
my $ny = int($nsamp/$nx + .5);
$nsamp = $nx * $ny;
# Compute scattering data using rtcontrib
my $cmd = "cnt $ndiv $ny $nx | rcalc -of -e '$tcal' " .
	"-e 'xp=(\$3+rand(.35*recno-15))*(($dim[1]-$dim[0])/$nx)+$dim[0]' " .
	"-e 'yp=(\$2+rand(.86*recno+11))*(($dim[3]-$dim[2])/$ny)+$dim[2]' " .
	"-e 'zp:$dim[4]-1e-5' " .
	q{-e 'Kbin=$1;x1=rand(1.21*recno+2.75);x2=rand(-3.55*recno-7.57)' } .
	q{-e '$1=xp;$2=yp;$3=zp;$4=Dx;$5=Dy;$6=Dz' } .
	"| rtcontrib -h -ff -n $nproc -c $nsamp -e '$kcal' -b kbin -bn $ndiv " .
	"-m $modnm -w -ab 5 -ad 700 -lw 3e-6 $octree " .
	"| rcalc -e '$tcal' " .
	"-e 'mod(n,d):n-floor(n/d)*d' -e 'Kbin=mod(recno-1,$ndiv)' " .
	q{-if3 -e '$1=(0.265*$1+0.670*$2+0.065*$3)/KprojOmega'};
my @darr = `$cmd`;
die "Failure running: $cmd\n" if ( $? );
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
printf "\t\t\t<Thickness unit=\"Meter\">%.3f</Thickness>\n", $dim[5] - $dim[4];
printf "\t\t\t<Width unit=\"Meter\">%.3f</Width>\n", $dim[1] - $dim[0];
printf "\t\t\t<Height unit=\"Meter\">%.3f</Height>\n", $dim[3] - $dim[2];
print "\t\t\t<DeviceType>Integral</DeviceType>\n";
# Output MGF description if requested
if ( $geout ) {
	print "\t\t\t<Geometry format=\"MGF\" unit=\"Meter\">\n";
	printf "xf -t %.6f %.6f 0\n", -($dim[0]+$dim[1])/2, -($dim[2]+$dim[3])/2;
	system "cat $mgfscn";
	print "xf\n";
	print "\t\t\t</Geometry>\n";
}
print '			</Material>
		<DataDefinition>
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
	<WavelengthData>
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
# Output computed data (transposed order)
for (my $od = 0; $od < $ndiv; $od++) {
	for (my $id = 0; $id < $ndiv; $id++) {
		print $darr[$ndiv*$id + $od];
	}
	print "\n";
}
# Output XML epilogue
print
'		</ScatteringData>
	</WavelengthDataBlock>
	</WavelengthData>
</Layer>
</Optical>
</WindowElement>
';
# Clean up temporary files
system "rm -rf $td";
