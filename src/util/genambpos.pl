#!/usr/bin/perl -w
# RCSid $Id: genambpos.pl,v 2.1 2014/04/24 23:15:42 greg Exp $
#
# Visualize ambient positions and gradients
#
use strict;
sub userror {
	print STDERR "Usage: genambpos [-l lvl][-s scale][-p][-d] scene.amb > ambloc.rad\n";
	exit 1;
}
my $lvlsel = -1;
my $scale = 0.25;
my $doposgrad = 0;
my $dodirgrad = 0;
my $savedARGV = "genambpos @ARGV";
# Get options
while ($#ARGV >= 0) {
	if ("$ARGV[0]" =~ /^-p/) {
		$doposgrad=1;
	} elsif ("$ARGV[0]" =~ /^-d/) {
		$dodirgrad=1;
	} elsif ("$ARGV[0]" =~ /^-l/) {
		$lvlsel = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-s/) {
		$scale = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-./) {
		userror();
	} else {
		last;
	}
	shift @ARGV;
}
userror() if ($#ARGV != 0);
my $cmd = "getinfo < $ARGV[0] " . 
		q[| sed -n 's/^.* -aa \([.0-9][.0-9]*\) .*$/\1/p'];
my $ambacc=`$cmd`;
die "Missing -aa setting in header\n" if (! $ambacc );
$scale *= $ambacc**.25;
my $outfmt = '
void glow posglow
0
0
4 ${agr} ${agg} ${agb} 0

posglow sphere position${recno}
0
0
4 ${px} ${py} ${pz} ${psiz}

posglow cone pgarrow${recno}
0
0
8
	${ cx0 }	${ cy0 }	${ cz0 }
	${ cx1 }	${ cy1 }	${ cz1 }
	${ cr0 }	0
';
my $posgradfmt = '
void brightfunc pgpat
2 posfunc ambpos.cal
0
6 ${ px } ${ py } ${ pz } ${ pgx } ${ pgy } ${ pgz }

pgpat glow pgval
0
0
4 ${avr} ${avg} ${avb} 0

void mixfunc pgeval
4 pgval void ellipstencil ambpos.cal
0
9 ${ px } ${ py } ${ pz } ${ux/r0} ${uy/r0} ${uz/r0} ${vx/r1} ${vy/r1} ${vz/r1}

pgeval polygon pgellipse${recno}
0
0
12
	${ px1 } ${ py1 } ${ pz1 }
	${ px2 } ${ py2 } ${ pz2 }
	${ px3 } ${ py3 } ${ pz3 }
	${ px4 } ${ py4 } ${ pz4 }
';
$outfmt .= $posgradfmt if ($doposgrad);
my $dirgradfmt='
void glow tipglow
0
0
4 ${2*agr} ${2*agg} ${2*agb} 0

tipglow sphere atip
0
0
4 ${ cx1 } ${ cy1 } ${ cz1 } ${psiz/7}

void brightfunc dgpat
2 dirfunc ambpos.cal
0
9 ${ px } ${ py } ${ pz } ${ nx } ${ ny } ${ nz } ${ dgx } ${ dgy } ${ dgz }

dgpat glow dgval
0
0
4 ${avr} ${avg} ${avb} 0

dgval ring dgdisk${recno}a
0
0
8
	${ px+dgx*.0001 } ${ py+dgy*.0001 } ${ pz+dgz*.0001 }
	${ dgx } ${ dgy } ${ dgz }
	0	${ r0/2 }

dgval ring dgdisk${recno}b
0
0
8
	${ px-dgx*.0001 } ${ py-dgy*.0001 } ${ pz-dgz*.0001 }
	${ -dgx } ${ -dgy } ${ -dgz }
	0	${ r0/2 }
';
$outfmt .= $dirgradfmt if ($dodirgrad);
# Load & convert ambient values
print "# Output produced by: $savedARGV\n";
system "lookamb -h -d $ARGV[0] | rcalc -e 'LV:$lvlsel;SF:$scale' -f rambpos.cal -o '$outfmt'\n";
exit;
