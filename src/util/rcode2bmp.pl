#!/usr/bin/perl
# RCSid $Id: rcode2bmp.pl,v 2.2 2020/01/21 18:00:17 greg Exp $
#
# Convert one or more rtpict outputs into BMP for convenient viewing
#

use strict;
use warnings;

my $xres=0;
my $yres=0;
my $pfilt="";

# Get options, which must precede input maps
while ($#ARGV >= 0 && "$ARGV[0]" =~ /^-[a-zA-Z]/) {
	if ("$ARGV[0]" eq '-x') {
		shift @ARGV;
		$xres = shift(@ARGV);
	} elsif ("$ARGV[0]" eq '-y') {
		shift @ARGV;
		$yres = shift(@ARGV);
	} else {
		die "Bad option: $ARGV[0]\n";
	}
}
# Set up resizing operation if requested
if ($xres > 0 && $yres > 0) {
	$pfilt = "pfilt -x $xres -y $yres -pa 1 -1 -r .6";
}
# Load/convert each file
INPUT_MAP:
while ($#ARGV >= 0) {
	if ("$ARGV[0]" =~ /\.bmp$/i) {
		print STDERR "Skipping file '$ARGV[0]' already in BMP format\n";
		shift @ARGV;
		next INPUT_MAP;
	}
	my $format = `getinfo < '$ARGV[0]' | sed -n 's/^FORMAT= *//p'`;
	chomp $format;
	die "Cannot get format from $ARGV[0]\n" if ( $? || ! $format );
	my ($dest) = ("$ARGV[0]" =~ /^([^.]+)/);
	$dest .= ".bmp";
	my $cmd="";
	if ("$format" =~ /^32-bit_rle_(rgb|xyz)e$/) {
		if ($pfilt) {
			$cmd = $pfilt . " '$ARGV[0]' | ra_bmp -e auto - '$dest'";
		} else {
			$cmd = "ra_bmp -e auto '$ARGV[0]' '$dest'";
		}
	} elsif ("$format" eq "16-bit_encoded_depth") {
		$cmd = "rcode_depth -r -ff -ho -Ho '$ARGV[0]' ";
		$cmd .= q{| rcalc -if -of -e 'cond=9e9-$1;$1=$1' | total -if -u};
		my $dmax=`$cmd`;
		$dmax = 2**(int(log($dmax)/log(2))+1);
		my $unit=`getinfo < '$ARGV[0]' | sed -n 's/^REFDEPTH= *[0-9.]*[^a-zA-Z]*//p'`;
		chomp $unit;
		$unit="Depth" if ( ! $unit );
		$cmd = "rcode_depth -r -ff '$ARGV[0]' | pvalue -r -df -b ";
		$cmd .= "| $pfilt " if ($pfilt);
		$cmd .= "| falsecolor -l '$unit' -m 1 -s $dmax | ra_bmp - '$dest'";
	} elsif ("$format" =~ /[1-9][0-9]*-bit_indexed_name$/) {
		$cmd = "rcode_ident -r -n '$ARGV[0]' " .
			"| getinfo +d -c rcalc -e 'cc(x):(.1+.8*rand(x))^2' " .
			q{-e '$1=cc(.398*$1-11.2);$2=cc(-1.152*$1+41.7);$3=cc(8.571*$1-8.15)' } .
			"| pvalue -r -d ";
		$cmd .= "| $pfilt " if ($pfilt);
		$cmd .= "| ra_bmp - '$dest'";
	} elsif ("$format" eq "32-bit_encoded_normal") {
		$cmd = "rcode_norm -r -ff '$ARGV[0]' | getinfo +d -c " .
			"rcalc -if3 -of -e `vwright v < '$ARGV[0]'` " .
				q{-e 'dot(vx,vy,vz)=vx*$1+vy*$2+vz*$3' } .
				"-e 'h=dot(vhx,vhy,vhz)' " .
				"-e 'v=dot(vvx,vvy,vvz)' " .
				"-e 'n=-dot(vdx,vdy,vdz)' " .
				"-f hsv_rgb.cal -e 'hue=mod(180/PI*atan2(v,h),360)' " .
				"-e 'sat=sqrt(h*h+v*v)' -e 'val=if(n,1,1+n)' " .
				q{-e '$1=R(hue,sat,val);$2=G(hue,sat,val);$3=B(hue,sat,val)' } .
				"| pvalue -r -df ";
		$cmd .= "| $pfilt " if ($pfilt);
		$cmd .= "| ra_bmp - '$dest'";
	}
	if ($cmd) {
		system $cmd;
		die "Cannot convert $ARGV[0]\n" if ( $? );
		print "Converted $ARGV[0] to $dest\n";
	} else {
		print STDERR "Skipping unsupported format '$format' for input '$ARGV[0]'\n";
	}
	shift @ARGV;
}
# EOF
