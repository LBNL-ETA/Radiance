#!/usr/bin/perl -w
# RCSid $Id: falsecolor.pl,v 2.1 2010/10/05 00:59:39 greg Exp $

use strict;
use File::Temp qw/ tempdir /;
use POSIX qw/ floor /;

my $mult = 179.0;    # Multiplier. Default W/sr/m2 -> cd/m2
my $label = 'cd/m2';    # Units shown in legend
my $scale = 1000;    # Top of the scale
my $decades = 0;    # Default is linear mapping
my $redv = 'def_red(v)';    # Mapping function for R,G,B
my $grnv = 'def_grn(v)';
my $bluv = 'def_blu(v)';
my $ndivs = 8;    # Number of lines in legend
my $picture = '-';
my $cpict = '';
my $legwidth = 100;   # Legend width and height
my $legheight = 200;
my $docont = '';    # Contours
my $loff = 0;    # Offset for drop-shadow
my $doextrem = 0;    # Don't mark extrema
my $needfile = 0;

while ($#ARGV >= 0) {
	# Options with qualifiers
    if ("$ARGV[0]" eq '-lw') {    # Legend width
        $legwidth = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-lh') {    # Legend height
        $legheight = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-m') {    # Multiplier
        $mult = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-s') {    # Scale
        $scale = $ARGV[1];
        shift @ARGV;
        if ($scale =~ m/[aA].*/) {
            $needfile = 1;
        }
    } elsif ("$ARGV[0]" eq '-l') {    # Label
        $label = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-log') {    # Logarithmic mapping
        $decades = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-r') {
        $redv = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-g') {
        $grnv = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-b') {
        $bluv = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-pal') {
	$redv = "$ARGV[1]_red(v)";
	$grnv = "$ARGV[1]_grn(v)";
	$bluv = "$ARGV[1]_blu(v)";
	shift @ARGV;
    } elsif ("$ARGV[0]" eq '-i') {    # Image for intensity mapping
        $picture = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-p') {    # Image for background
        $cpict = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-ip' || "$ARGV[0]" eq '-pi') {
        $picture = $ARGV[1];
        $cpict = $ARGV[1];
        shift @ARGV;
    } elsif ("$ARGV[0]" eq '-n') {    # Number of contour lines
        $ndivs = $ARGV[1];
        shift @ARGV;

	# Switches
    } elsif ("$ARGV[0]" eq '-cl') {    # Contour lines
        $docont = 'a';
        $loff = 12;
    } elsif ("$ARGV[0]" eq '-cb') {    # Contour bands
        $docont = 'b';
        $loff = 13;
    } elsif ("$ARGV[0]" eq '-e') {
        $doextrem = 1;
        $needfile = 1;

	# Oops! Illegal option
    } else {
        die("bad option \"$ARGV[0]\"\n");
    }
    shift @ARGV;
}

# Temporary directory. Will be removed upon successful program exit.
my $td = tempdir( CLEANUP => 1 );

if ($needfile == 1 && $picture eq '-') {
    $picture = $td . '/' . $picture;
    open(FHpic, ">$picture") or
            die("Unable to write to file $picture\n");
    close(FHpic);
}

# Find a good scale for auto mode.
if ($scale =~ m/[aA].*/) {
    my $LogLmax = `phisto $picture| tail -2| sed -n '1s/[0-9]*\$//p'`;
    $scale = $mult / 179 * 10**$LogLmax;
}

my $pc0 = $td . '/pc0.cal';
open(FHpc0, ">$pc0");
print FHpc0 <<EndOfPC0;
PI : 3.14159265358979323846 ;
scale : $scale ;
mult : $mult ;
ndivs : $ndivs ;
gamma : 2.2;

or(a,b) : if(a,a,b);
EPS : 1e-7;
neq(a,b) : if(a-b-EPS,1,b-a-EPS);
btwn(a,x,b) : if(a-x,-1,b-x);
clip(x) : if(x-1,1,if(x,x,0));
frac(x) : x - floor(x);
boundary(a,b) : neq(floor(ndivs*a+.5),floor(ndivs*b+.5));

spec_red(x) = 1.6*x - .6;
spec_grn(x) = if(x-.375, 1.6-1.6*x, 8/3*x);
spec_blu(x) = 1 - 8/3*x;

pm3d_red(x) = sqrt(x) ^ gamma;
pm3d_grn(x) = x*x*x ^ gamma;
pm3d_blu(x) = clip(sin(2*PI*x)) ^ gamma;

hot_red(x) = clip(3*x) ^ gamma;
hot_grn(x) = clip(3*x - 1) ^ gamma;
hot_blu(x) = clip(3*x - 2) ^ gamma;

interp_arr2(i,x,f):(i+1-x)*f(i)+(x-i)*f(i+1);
interp_arr(x,f):if(x-1,if(f(0)-x,interp_arr2(floor(x),x,f),f(f(0))),f(1));

def_redp(i):select(i,0.18848,0.05468174,
0.00103547,8.311144e-08,7.449763e-06,0.0004390987,0.001367254,
0.003076,0.01376382,0.06170773,0.1739422,0.2881156,0.3299725,
0.3552663,0.372552,0.3921184,0.4363976,0.6102754,0.7757267,
0.9087369,1,1,0.9863);
def_red(x):interp_arr(x/0.0454545+1,def_redp);

def_grnp(i):select(i,0.0009766,2.35501e-05,
0.0008966244,0.0264977,0.1256843,0.2865799,0.4247083,0.4739468,
0.4402732,0.3671876,0.2629843,0.1725325,0.1206819,0.07316644,
0.03761026,0.01612362,0.004773749,6.830967e-06,0.00803605,
0.1008085,0.3106831,0.6447838,0.9707);
def_grn(x):interp_arr(x/0.0454545+1,def_grnp);

def_blup(i):select(i,0.2666,0.3638662,0.4770437,
0.5131397,0.5363797,0.5193677,0.4085123,0.1702815,0.05314236,
0.05194055,0.08564082,0.09881395,0.08324373,0.06072902,
0.0391076,0.02315354,0.01284458,0.005184709,0.001691774,
2.432735e-05,1.212949e-05,0.006659406,0.02539);
def_blu(x):interp_arr(x/0.0454545+1,def_blup);

isconta = if(btwn(0,v,1),or(boundary(vleft,vright),boundary(vabove,vbelow)),-1);
iscontb = if(btwn(0,v,1),btwn(.4,frac(ndivs*v),.6),-1);

ra = 0;
ga = 0;
ba = 0;

in = 1;

ro = if(in,clip($redv),ra);
go = if(in,clip($grnv),ga);
bo = if(in,clip($bluv),ba);
EndOfPC0

my $pc1 = $td . '/pc1.cal';
open(FHpc1, ">$pc1");
print FHpc1 <<EndOfPC1;
norm : mult/scale/le(1);

v = map(li(1)*norm);

vleft = map(li(1,-1,0)*norm);
vright = map(li(1,1,0)*norm);
vabove = map(li(1,0,1)*norm);
vbelow = map(li(1,0,-1)*norm);

map(x) = x;

ra = ri(nfiles);
ga = gi(nfiles);
ba = bi(nfiles);
EndOfPC1

my $pc0args = "-f $pc0";
my $pc1args = "-f $pc1";

# Contour lines or bands
if ($docont ne '') {
    $pc0args .= " -e 'in=iscont$docont'";
}

if ($cpict eq '') {
    $pc1args .= " -e 'ra=0;ga=0;ba=0'";
} elsif ($cpict eq $picture) {
    $cpict = '';
}

# Logarithmic mapping
if ($decades > 0) {
    $pc1args .= " -e 'map(x)=if(x-10^-$decades,log10(x)/$decades+1,0)'";
}

# Colours in the legend
my $scolpic = $td . '/scol.hdr';
# Labels in the legend
my $slabpic = $td . '/slab.hdr';
my $cmd;

if ($legwidth > 20 && $legheight > 40) {
	# Legend: Create the background colours
    $cmd = "pcomb $pc0args -e 'v=(y+.5)/yres;vleft=v;vright=v'";
    $cmd .= " -e 'vbelow=(y-.5)/yres;vabove=(y+1.5)/yres'";
    $cmd .= " -x $legwidth -y $legheight > $scolpic";
    system $cmd;

	# Legend: Create the text labels
    my $sheight = floor($legheight / $ndivs + 0.5);
    my $text = "$label";
    my $value;
    my $i;
    for ($i=0; $i<$ndivs; $i++) {
        my $imap = ($ndivs - 0.5 - $i) / $ndivs;
        if ($decades > 0) {
            $value = $scale * 10**(($imap - 1) * $decades);
        } else {
            $value = $scale * $imap;
        }
        # Have no more than 3 decimal places
        $value =~ s/(\.[0-9]{3})[0-9]*/$1/;
        $text .= "\n$value";
    }
    $cmd = "echo '$text' | psign -s -.15 -cf 1 1 1 -cb 0 0 0";
    $cmd .= " -h $sheight > $slabpic";
    system $cmd;
} else {
	# Legend is too small to be legible. Don't bother doing one.
    $legwidth = 0;
    $legheight = 0;
    open(FHscolpic, ">$scolpic");
    print FHscolpic '\n-Y 1 +X 1\naaa\n';
    close(FHscolpic);
    open(FHslabpic, ">$slabpic");
    print FHslabpic '\n-Y 1 +X 1\naaa\n';
    close(FHslabpic);
}

my $loff1 = $loff - 1;
# Command line without extrema
$cmd = "pcomb $pc0args $pc1args $picture $cpict";
$cmd .= "| pcompos $scolpic 0 0 +t .1";
$cmd .= " \"\!pcomb -e 'lo=1-gi(1)' $slabpic\" 2 $loff1";
$cmd .= " -t .5 $slabpic 0 $loff - $legwidth 0";

if ($doextrem == 1) {
	# Get min/max image luminance
    my $cmd1 = 'pextrem -o ' . $picture;
    my $retval = `$cmd1`;
    # e.g.
    # x   y   r            g            b
    # 193 207 3.070068e-02 3.118896e-02 1.995850e-02
    # 211 202 1.292969e+00 1.308594e+00 1.300781e+00

    my @extrema = split(/\s/, $retval);
    my $lxmin = $extrema[0] + $legwidth;
    my $ymin = $extrema[1];
    my $rmin = $extrema[2];
    my $gmin = $extrema[3];
    my $bmin = $extrema[4];
    my $lxmax = $extrema[5] + $legwidth;
    my $ymax = $extrema[6];
    my $rmax = $extrema[7];
    my $gmax = $extrema[8];
    my $bmax = $extrema[9];

	# Weighted average of R,G,B
    my $minpos = "$lxmin $ymin";
    my $minval = ($rmin * .27 + $gmin * .67 + $bmin * .06) * $mult;
    $minval =~ s/(\.[0-9]{3})[0-9]*/$1/;
    my $maxpos = "$lxmax $ymax";
    my $maxval = ($rmax * .27 + $gmax * .67 + $bmax * .06) * $mult;
    $maxval =~ s/(\.[0-9]{3})[0-9]*/$1/;

	# Create the labels for min/max intensity
    my $minvpic = $td . '/minv.hdr';
    $cmd1 = "psign -s -.15 -a 2 -h 16 $minval > $minvpic";
    system $cmd1;
    my $maxvpic = $td . '/maxv.hdr';
    $cmd1 = "psign -s -.15 -a 2 -h 16 $maxval > $maxvpic";
    system $cmd1;

	# Add extrema labels to command line
    $cmd .= " $minvpic $minpos $maxvpic $maxpos";
}

# Process image and combine with legend
system "$cmd";

#EOF
