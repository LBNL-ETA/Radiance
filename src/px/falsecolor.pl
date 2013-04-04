#!/usr/bin/perl -w
# RCSid $Id: falsecolor.pl,v 2.9 2013/04/04 02:59:20 greg Exp $

use warnings;
use strict;
use File::Temp qw/ tempdir /;
use POSIX qw/ floor /;

my @palettes = ('def', 'spec', 'pm3d', 'hot', 'eco');

my $mult = 179.0;              # Multiplier. Default W/sr/m2 -> cd/m2
my $label = 'cd/m2';           # Units shown in legend
my $scale = 1000;              # Top of the scale
my $decades = 0;               # Default is linear mapping
my $pal = 'def';               # Palette
my $redv = "${pal}_red(v)";    # Mapping functions for R,G,B
my $grnv = "${pal}_grn(v)";
my $bluv = "${pal}_blu(v)";
my $ndivs = 8;                 # Number of lines in legend
my $picture = '-';
my $cpict = '';
my $legwidth = 100;            # Legend width and height
my $legheight = 200;
my $docont = '';               # Contours: -cl and -cb
my $doposter = 0;              # Posterization: -cp
my $loff = 0;                  # Offset to align with values
my $doextrem = 0;              # Don't mark extrema
my $needfile = 0;
my $showpal = 0;               # Show availabel colour palettes

while ($#ARGV >= 0) {
    $_ = shift;
    # Options with qualifiers
    if (m/-lw/) {              # Legend width
        $legwidth = shift;
    } elsif (m/-lh/) {         # Legend height
        $legheight = shift;
    } elsif (m/-m/) {          # Multiplier
        $mult = shift;
    } elsif (m/-spec/) {
        die("depricated option '-spec'. Please use '-pal spec' instead.");
    } elsif (m/-s/) {          # Scale
        $scale = shift;
        if ($scale =~ m/[aA].*/) {
            $needfile = 1;
        }
    } elsif (m/-l$/) {         # Label
        $label = shift;
    } elsif (m/-log/) {        # Logarithmic mapping
        $decades = shift;
    } elsif (m/-r/) {          # Custom palette functions for R,G,B
        $redv = shift;
    } elsif (m/-g/) {
        $grnv = shift;
    } elsif (m/-b/) {
        $bluv = shift;
    } elsif (m/-pal$/) {        # Color palette
        $pal = shift;
        if (! grep $_ eq $pal, @palettes) {
            die("invalid palette '$pal'.\n");
        }
        $redv = "${pal}_red(v)";
        $grnv = "${pal}_grn(v)";
        $bluv = "${pal}_blu(v)";
    } elsif (m/-i$/) {          # Image for intensity mapping
        $picture = shift;
    } elsif (m/-p$/) {         # Image for background
        $cpict = shift;
    } elsif (m/-ip/ || m/-pi/) {
        $picture = shift;
        $cpict = $picture;
    } elsif (m/-n/) {          # Number of contour lines
        $ndivs = shift;

    # Switches
    } elsif (m/-cl/) {         # Contour lines
        $docont = 'a';
        $loff = 0.48;
    } elsif (m/-cb/) {         # Contour bands
        $docont = 'b';
        $loff = 0.52;
    } elsif (m/-cp/) {              # Posterize
        $doposter = 1;
    } elsif (m/-palettes/) {        # Show all available palettes
        $scale   = 45824;           # 256 * 179
        $showpal = 1;
    } elsif (m/-e/) {
        $doextrem = 1;
        $needfile = 1;

    # Oops! Illegal option
    } else {
        die("bad option \"$_\"\n");
    }
}

# Temporary directory. Will be removed upon successful program exit.
my $td = tempdir( CLEANUP => 1 );

if ($needfile == 1 && $picture eq '-') {
    # Pretend that $td/stdin.rad is the actual filename.
    $picture = "$td/stdin.hdr";
    open(FHpic, ">$picture") or
            die("Unable to write to file $picture\n");
    # Input is from STDIN: Capture to file.
    while (<>) {
        print FHpic;
    }
    close(FHpic);

    if ($cpict eq '-') {
        $cpict = "$td/stdin.hdr";
    }
}

# Find a good scale for auto mode.
if ($scale =~ m/[aA].*/) {
    my @histo = split(/\s/, `phisto $picture| tail -2`);
    # e.g. $ phisto tests/richmond.hdr| tail -2
    # 3.91267	14
    # 3.94282	6
    my $LogLmax = $histo[0];
    $scale = $mult / 179 * 10**$LogLmax;
}

my $pc0 = "$td/pc0.cal";
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
pm3d_grn(x) = (x*x*x) ^ gamma;
pm3d_blu(x) = clip(sin(2*PI*clip(x))) ^ gamma;

hot_red(x) = clip(3*x) ^ gamma;
hot_grn(x) = clip(3*x - 1) ^ gamma;
hot_blu(x) = clip(3*x - 2) ^ gamma;

eco_red(x) = clip(2*x) ^ gamma;
eco_grn(x) = clip(2*(x-0.5)) ^ gamma;
eco_blu(x) = clip(2*(0.5-x)) ^ gamma;

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
close FHpc0;

my $pc1 = "$td/pc1.cal";
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
close FHpc1;

my $pc0args = "-f $pc0";
my $pc1args = "-f $pc1";

if ($showpal == 1) {
    my $pc = "pcompos -a 1";
    foreach my $pal (@palettes) {
        my $fcimg = "$td/$pal.hdr";
        my $lbimg = "$td/${pal}_label.hdr";
        system "psign -cb 0 0 0 -cf 1 1 1 -h 20 $pal > $lbimg";

        my $cmd = qq[pcomb $pc0args -e "v=x/256"];
        $cmd .= qq[ -e "ro=clip(${pal}_red(v));go=clip(${pal}_grn(v));bo=clip(${pal}_blu(v))"];
        $cmd .= qq[ -x 256 -y 30 > $fcimg];
        system "$cmd";
        $pc .= " $fcimg $lbimg";
    }
    system "$pc";
    exit 0;
}

# Contours
if ($docont ne '') {
    # -cl -> $docont = a
    # -cb -> $docont = b
    $pc0args .= qq[ -e "in=iscont$docont"];
} elsif ($doposter == 1) {
    # -cp -> $doposter = 1
    $pc0args .= qq[ -e "ro=${pal}_red(seg(v));go=${pal}_grn(seg(v));bo=${pal}_blu(seg(v))"];
    $pc0args .= q[ -e "seg(x)=(floor(v*ndivs)+.5)/ndivs"];
}

if ($cpict eq '') {
    $pc1args .= qq[ -e "ra=0;ga=0;ba=0"];
} elsif ($cpict eq $picture) {
    $cpict = '';
}

# Logarithmic mapping
if ($decades > 0) {
    $pc1args .= qq[ -e "map(x)=if(x-10^-$decades,log10(x)/$decades+1,0)"];
}

# Colours in the legend
my $scolpic = "$td/scol.hdr";

# Labels in the legend
my $slabpic = "$td/slab.hdr";
my $cmd;

if (($legwidth > 20) && ($legheight > 40)) {
    # Legend: Create the text labels
    my $sheight = floor($legheight / $ndivs + 0.5);
    $legheight = $sheight * $ndivs;
    $loff = floor($loff * $sheight + 0.5);
    my $text = "$label";
    for (my $i=0; $i<$ndivs; $i++) {
        my $imap = ($ndivs - 0.5 - $i) / $ndivs;
        my $value = $scale;
        if ($decades > 0) {
            $value *= 10**(($imap - 1) * $decades);
        } else {
            $value *= $imap;
        }

        # Have no more than 3 decimal places
        $value =~ s/(\.[0-9]{3})[0-9]*/$1/;
        $text .= "\n$value";
    }
    open PSIGN, "| psign -s -.15 -cf 1 1 1 -cb 0 0 0 -h $sheight > $slabpic";
    print PSIGN "$text\n";
    close PSIGN;

    # Legend: Create the background colours
    $cmd = qq[pcomb $pc0args];
    $cmd .= q[ -e "v=(y+.5)/yres;vleft=v;vright=v"];
    $cmd .= q[ -e "vbelow=(y-.5)/yres;vabove=(y+1.5)/yres"];
    $cmd .= qq[ -x $legwidth -y $legheight > $scolpic];
    system $cmd;
} else {
    # Legend is too small to be legible. Don't bother doing one.
    $legwidth = 0;
    $legheight = 0;
    $loff = 0;

    # Create dummy colour scale and legend labels so we don't
    # need to change the final command line.
    open(FHscolpic, ">$scolpic");
    print FHscolpic "\n-Y 1 +X 1\naaa\n";
    close(FHscolpic);
    open(FHslabpic, ">$slabpic");
    print FHslabpic "\n-Y 1 +X 1\naaa\n";
    close(FHslabpic);
}

# Legend: Invert the text labels (for dropshadow)
my $slabinvpic = "$td/slabinv.hdr";
$cmd = qq[pcomb -e "lo=1-gi(1)" $slabpic > $slabinvpic];
system $cmd;

my $loff1 = $loff - 1;

# Command line without extrema
$cmd = qq[pcomb $pc0args $pc1args $picture $cpict];
$cmd .= qq[ | pcompos $scolpic 0 0 +t .1 $slabinvpic 2 $loff1];
$cmd .= qq[ -t .5 $slabpic 0 $loff - $legwidth 0];

if ($doextrem == 1) {
    # Get min/max image luminance
    my $cmd1 = 'pextrem -o ' . $picture;
    my $retval = `$cmd1`;
    # e.g.
    # x   y   r            g            b
    # 193 207 3.070068e-02 3.118896e-02 1.995850e-02
    # 211 202 1.292969e+00 1.308594e+00 1.300781e+00

    my @extrema = split(/\s/, $retval);
    my ($lxmin, $ymin, $rmin, $gmin, $bmin, $lxmax, $ymax, $rmax, $gmax, $bmax) = @extrema;
    $lxmin += $legwidth;
    $lxmax += $legwidth;

    # Weighted average of R,G,B
    my $minpos = "$lxmin $ymin";
    my $minval = ($rmin * .27 + $gmin * .67 + $bmin * .06) * $mult;
    $minval =~ s/(\.[0-9]{3})[0-9]*/$1/;
    my $maxval = ($rmax * .27 + $gmax * .67 + $bmax * .06) * $mult;
    $maxval =~ s/(\.[0-9]{3})[0-9]*/$1/;

    # Create the labels for min/max intensity
    my $minvpic = "$td/minv.hdr";
    system "psign -s -.15 -a 2 -h 16 $minval > $minvpic";
    my $maxvpic = "$td/maxv.hdr";
    system "psign -s -.15 -a 2 -h 16 $maxval > $maxvpic";

    # Add extrema labels to command line
    $cmd .= qq[ $minvpic $minpos $maxvpic $lxmax $ymax];
}

# Process image and combine with legend
system "$cmd";

#EOF
