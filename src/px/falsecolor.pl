#!/usr/bin/perl -w
# RCSid $Id: falsecolor.pl,v 2.14 2019/08/25 22:24:01 greg Exp $

use warnings;
use strict;
use File::Temp qw/ tempdir /;
use POSIX qw/ floor /;

my @palettes = ('def', 'spec', 'pm3d', 'hot', 'eco', 'tbo');

my $mult = 179.0;              # Multiplier. Default W/sr/m2 -> cd/m2
my $label = 'cd/m2';           # Units shown in legend
my $scale = 1000;              # Top of the scale
my $scaledigits = 3;           # Number of maximum digits for numbers in legend
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
my $haszero = 1;               # print 0 in scale for falsecolor images only
my $docont = '';               # Contours: -cl and -cb
my $doposter = 0;              # Posterization: -cp
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
    } elsif (m/-d/) {          # Number of maximum digits for numbers in legend
        $scaledigits = shift;
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
        $haszero = 0;
    } elsif (m/-cb/) {         # Contour bands
        $docont = 'b';
        $haszero = 0;
    } elsif (m/-cp/) {              # Posterize
        $doposter = 1;
        $haszero = 0;
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

if (($legwidth <= 20) || ($legheight <= 40)) {
    # Legend is too small to be legible. Don't bother doing one.
    $legwidth = 0;
    $legheight = 0;
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

if ($docont ne '') {
    # -cl -> $docont = a
    # -cb -> $docont = b
    my $newv = join( "(v-1/ndivs)*ndivs/(ndivs-1)", split("v", $redv) );
    $redv = $newv;
    $newv = join( "(v-1/ndivs)*ndivs/(ndivs-1)", split("v", $bluv) );
    $bluv = $newv;
    $newv = join( "(v-1/ndivs)*ndivs/(ndivs-1)", split("v", $grnv) );
    $grnv = $newv;
} elsif ($doposter == 1) {
    # -cp -> $doposter = 1
    my $newv = join( "seg2(v)", split("v", $redv) );
    $redv = $newv;
    $newv = join( "seg2(v)", split("v", $bluv) );
    $bluv = $newv;
    $newv = join( "seg2(v)", split("v", $grnv) );
    $grnv = $newv;
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
boundary(a,b) : neq(floor(ndivs*a),floor(ndivs*b));
round(x):if(x-floor(x)-0.5,ceil(x),floor(x));

spec_red(x) = 1.6*x - .6;
spec_grn(x) = if(x-.375, 1.6-1.6*x, 8/3*x);
spec_blu(x) = 1 - 8/3*x;

pm3d_red(x) = sqrt(clip(x)) ^ gamma;
pm3d_grn(x) = clip(x*x*x) ^ gamma;
pm3d_blu(x) = clip(sin(2*PI*clip(x))) ^ gamma;

hot_red(x) = clip(3*x) ^ gamma;
hot_grn(x) = clip(3*x - 1) ^ gamma;
hot_blu(x) = clip(3*x - 2) ^ gamma;

eco_red(x) = clip(2*x) ^ gamma;
eco_grn(x) = clip(2*(x-0.5)) ^ gamma;
eco_blu(x) = clip(2*(0.5-x)) ^ gamma;

interp_arr2(i,x,f):(i+1-x)*f(i)+(x-i)*f(i+1);
interp_arr(x,f):if(x-1,if(f(0)-x,interp_arr2(floor(x),x,f),f(f(0))),f(1));

interp_tbo2(x,f):f(round(x*255)+1)+(f(min(255,round(x*255)+1)+1)-f(round(x*255)+1))*(x*255-round(x*255));
interp_tbo(x,f):if(x-1,f(256),if(x,interp_tbo2(x,f),f(1)));

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

tbo_redp(i):select(i,0.18995,0.19483,0.19956,0.20415,0.2086,0.21291,0.21708,0.22111,0.225,0.22875,0.23236,0.23582,0.23915,0.24234,0.24539,0.2483,0.25107,0.25369,0.25618,0.25853,0.26074,0.2628,0.26473,0.26652,0.26816,0.26967,0.27103,0.27226,0.27334,0.27429,0.27509,0.27576,0.27628,0.27667,0.27691,0.27701,0.27698,0.2768,0.27648,0.27603,0.27543,0.27469,0.27381,0.27273,0.27106,0.26878,0.26592,0.26252,0.25862,0.25425,0.24946,0.24427,0.23874,0.23288,0.22676,0.22039,0.21382,0.20708,0.20021,0.19326,0.18625,0.17923,0.17223,0.16529,0.15844,0.15173,0.14519,0.13886,0.13278,0.12698,0.12151,0.11639,0.11167,0.10738,0.10357,0.10026,0.0975,0.09532,0.09377,0.09287,0.09267,0.0932,0.09451,0.09662,0.09958,0.10342,0.10815,0.11374,0.12014,0.12733,0.13526,0.14391,0.15323,0.16319,0.17377,0.18491,0.19659,0.20877,0.22142,0.23449,0.24797,0.2618,0.27597,0.29042,0.30513,0.32006,0.33517,0.35043,0.36581,0.38127,0.39678,0.41229,0.42778,0.44321,0.45854,0.47375,0.48879,0.50362,0.51822,0.53255,0.54658,0.56026,0.57357,0.58646,0.59891,0.61088,0.62233,0.63323,0.64362,0.65394,0.66428,0.67462,0.68494,0.69525,0.70553,0.71577,0.72596,0.7361,0.74617,0.75617,0.76608,0.77591,0.78563,0.79524,0.80473,0.8141,0.82333,0.83241,0.84133,0.8501,0.85868,0.86709,0.8753,0.88331,0.89112,0.8987,0.90605,0.91317,0.92004,0.92666,0.93301,0.93909,0.94489,0.95039,0.9556,0.96049,0.96507,0.96931,0.97323,0.97679,0.98,0.98289,0.98549,0.98781,0.98986,0.99163,0.99314,0.99438,0.99535,0.99607,0.99654,0.99675,0.99672,0.99644,0.99593,0.99517,0.99419,0.99297,0.99153,0.98987,0.98799,0.9859,0.9836,0.98108,0.97837,0.97545,0.97234,0.96904,0.96555,0.96187,0.95801,0.95398,0.94977,0.94538,0.94084,0.93612,0.93125,0.92623,0.92105,0.91572,0.91024,0.90463,0.89888,0.89298,0.88691,0.88066,0.87422,0.8676,0.86079,0.8538,0.84662,0.83926,0.83172,0.82399,0.81608,0.80799,0.79971,0.79125,0.7826,0.77377,0.76476,0.75556,0.74617,0.73661,0.72686,0.71692,0.7068,0.6965,0.68602,0.67535,0.66449,0.65345,0.64223,0.63082,0.61923,0.60746,0.5955,0.58336,0.57103,0.55852,0.54583,0.53295,0.51989,0.50664,0.49321,0.4796);
tbo_red(x):interp_tbo(x,tbo_redp) ^ gamma;

tbo_grnp(i):select(i,0.07176,0.08339,0.09498,0.10652,0.11802,0.12947,0.14087,0.15223,0.16354,0.17481,0.18603,0.1972,0.20833,0.21941,0.23044,0.24143,0.25237,0.26327,0.27412,0.28492,0.29568,0.30639,0.31706,0.32768,0.33825,0.34878,0.35926,0.3697,0.38008,0.39043,0.40072,0.41097,0.42118,0.43134,0.44145,0.45152,0.46153,0.47151,0.48144,0.49132,0.50115,0.51094,0.52069,0.5304,0.54015,0.54995,0.55979,0.56967,0.57958,0.5895,0.59943,0.60937,0.61931,0.62923,0.63913,0.64901,0.65886,0.66866,0.67842,0.68812,0.69775,0.70732,0.7168,0.7262,0.73551,0.74472,0.75381,0.76279,0.77165,0.78037,0.78896,0.7974,0.80569,0.81381,0.82177,0.82955,0.83714,0.84455,0.85175,0.85875,0.86554,0.87211,0.87844,0.88454,0.8904,0.896,0.90142,0.90673,0.91193,0.91701,0.92197,0.9268,0.93151,0.93609,0.94053,0.94484,0.94901,0.95304,0.95692,0.96065,0.96423,0.96765,0.97092,0.97403,0.97697,0.97974,0.98234,0.98477,0.98702,0.98909,0.99098,0.99268,0.99419,0.99551,0.99663,0.99755,0.99828,0.99879,0.9991,0.99919,0.99907,0.99873,0.99817,0.99739,0.99638,0.99514,0.99366,0.99195,0.98999,0.98775,0.98524,0.98246,0.97941,0.9761,0.97255,0.96875,0.9647,0.96043,0.95593,0.95121,0.94627,0.94113,0.93579,0.93025,0.92452,0.91861,0.91253,0.90627,0.89986,0.89328,0.88655,0.87968,0.87267,0.86553,0.85826,0.85087,0.84337,0.83576,0.82806,0.82025,0.81236,0.80439,0.79634,0.78823,0.78005,0.77181,0.76352,0.75519,0.74682,0.73842,0.73,0.7214,0.7125,0.7033,0.69382,0.68408,0.67408,0.66386,0.65341,0.64277,0.63193,0.62093,0.60977,0.59846,0.58703,0.57549,0.56386,0.55214,0.54036,0.52854,0.51667,0.50479,0.49291,0.48104,0.4692,0.4574,0.44565,0.43399,0.42241,0.41093,0.39958,0.38836,0.37729,0.36638,0.35566,0.34513,0.33482,0.32473,0.31489,0.3053,0.29599,0.28696,0.27824,0.26981,0.26152,0.25334,0.24526,0.2373,0.22945,0.2217,0.21407,0.20654,0.19912,0.19182,0.18462,0.17753,0.17055,0.16368,0.15693,0.15028,0.14374,0.13731,0.13098,0.12477,0.11867,0.11268,0.1068,0.10102,0.09536,0.0898,0.08436,0.07902,0.0738,0.06868,0.06367,0.05878,0.05399,0.04931,0.04474,0.04028,0.03593,0.03169,0.02756,0.02354,0.01963,0.01583);
tbo_grn(x):interp_tbo(x,tbo_grnp) ^ gamma;

tbo_blup(i):select(i,0.23217,0.26149,0.29024,0.31844,0.34607,0.37314,0.39964,0.42558,0.45096,0.47578,0.50004,0.52373,0.54686,0.56942,0.59142,0.61286,0.63374,0.65406,0.67381,0.693,0.71162,0.72968,0.74718,0.76412,0.7805,0.79631,0.81156,0.82624,0.84037,0.85393,0.86692,0.87936,0.89123,0.90254,0.91328,0.92347,0.93309,0.94214,0.95064,0.95857,0.96594,0.97275,0.97899,0.98461,0.9893,0.99303,0.99583,0.99773,0.99876,0.99896,0.99835,0.99697,0.99485,0.99202,0.98851,0.98436,0.97959,0.97423,0.96833,0.9619,0.95498,0.94761,0.93981,0.93161,0.92305,0.91416,0.90496,0.8955,0.8858,0.8759,0.86581,0.85559,0.84525,0.83484,0.82437,0.81389,0.80342,0.79299,0.78264,0.7724,0.7623,0.75237,0.74265,0.73316,0.72393,0.715,0.70599,0.69651,0.6866,0.67627,0.66556,0.65448,0.64308,0.63137,0.61938,0.60713,0.59466,0.58199,0.56914,0.55614,0.54303,0.52981,0.51653,0.50321,0.48987,0.47654,0.46325,0.45002,0.43688,0.42386,0.41098,0.39826,0.38575,0.37345,0.3614,0.34963,0.33816,0.32701,0.31622,0.30581,0.29581,0.28623,0.27712,0.26849,0.26038,0.2528,0.24579,0.23937,0.23356,0.22835,0.2237,0.2196,0.21602,0.21294,0.21032,0.20815,0.2064,0.20504,0.20406,0.20343,0.20311,0.2031,0.20336,0.20386,0.20459,0.20552,0.20663,0.20788,0.20926,0.21074,0.2123,0.21391,0.21555,0.21719,0.2188,0.22038,0.22188,0.22328,0.22456,0.2257,0.22667,0.22744,0.228,0.22831,0.22836,0.22811,0.22754,0.22663,0.22536,0.22369,0.22161,0.21918,0.2165,0.21358,0.21043,0.20706,0.20348,0.19971,0.19577,0.19165,0.18738,0.18297,0.17842,0.17376,0.16899,0.16412,0.15918,0.15417,0.1491,0.14398,0.13883,0.13367,0.12849,0.12332,0.11817,0.11305,0.10797,0.10294,0.09798,0.0931,0.08831,0.08362,0.07905,0.07461,0.07031,0.06616,0.06218,0.05837,0.05475,0.05134,0.04814,0.04516,0.04243,0.03993,0.03753,0.03521,0.03297,0.03082,0.02875,0.02677,0.02487,0.02305,0.02131,0.01966,0.01809,0.0166,0.0152,0.01387,0.01264,0.01148,0.01041,0.00942,0.00851,0.00769,0.00695,0.00629,0.00571,0.00522,0.00481,0.00449,0.00424,0.00408,0.00401,0.00401,0.0041,0.00427,0.00453,0.00486,0.00529,0.00579,0.00638,0.00705,0.0078,0.00863,0.00955,0.01055);
tbo_blu(x):interp_tbo(x,tbo_blup) ^ gamma;

isconta = if(btwn(1/ndivs/2,v,1+1/ndivs/2),or(boundary(vleft,vright),boundary(vabove,vbelow)),-1);
iscontb = if(btwn(1/ndivs/2,v,1+1/ndivs/2),-btwn(.1,frac(ndivs*v),.9),-1);

seg(x)=(floor(v*ndivs)+.5)/ndivs;
seg2(x)=(seg(x)-1/ndivs)*ndivs/(ndivs-1);

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
if ($legwidth > 0) {
    # Legend: Create the text labels
    my $sheight = floor($legheight / $ndivs );
    my $theight = floor($legwidth/(8/1.67));
    my $stheight = $sheight <= $theight ? $sheight : $theight;
    my $vlegheight = $sheight * $ndivs * (1+1.5/$ndivs);
    my $text = "$label\n";
    my $tslabpic = "$td/slabT.hdr";
    open PSIGN, "| psign -s -.15 -cf 1 1 1 -cb 0 0 0 -h $stheight > $tslabpic";
    print PSIGN "$text";
    close PSIGN;
    my $loop = $ndivs+$haszero;
    my $hlegheight = $sheight * ($loop) + $sheight * .5;
    my $pcompost = qq[pcompos -b 0 0 0 =-0 $tslabpic 0 $hlegheight ];
    for (my $i=0; $i<$loop; $i++) {
        my $imap = ($ndivs - $i) / $ndivs;
        my $value = $scale;
        if ($decades > 0) {
            $value *= 10**(($imap - 1) * $decades);
        } else {
            $value *= $imap;
        }
        # Have no more than 3 decimal places
        $value =~ s/(\.[0-9]{$scaledigits})[0-9]*/$1/;
        $text = "$value\n";
        $tslabpic = "$td/slab$i.hdr";
        open PSIGN, "| psign -s -.15 -cf 1 1 1 -cb 0 0 0 -h $stheight > $tslabpic";
        print PSIGN "$text";
        close PSIGN;
        $hlegheight = $sheight * ($loop - $i - 1) + $sheight * .5;
        $pcompost .= qq[=-0 $tslabpic 0 $hlegheight ];
    }
    $pcompost .= qq[ > $slabpic];
    system $pcompost;

    # Legend: Create the background colours
    $cmd = qq[pcomb $pc0args];
    $cmd .= qq[ -x $legwidth -y $vlegheight];
    $cmd .= qq[ -e "v=(y+.5-$sheight)/(yres/(1+1.5/$ndivs));vleft=v;vright=v"];
    $cmd .= qq[ -e "vbelow=(y-.5-$sheight)/(yres/(1+1.5/$ndivs));vabove=(y+1.5-$sheight)/(yres/(1+1.5/$ndivs))"];
    $cmd .= qq[ -e "ra=0;ga=0;ba=0;"];
    $cmd .= qq[ > $scolpic];
    system $cmd;
} else {
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


my $sh0 = -floor($legheight / $ndivs / 2);
if ($haszero < 1) {
    $sh0 = -floor($legheight / ($ndivs)*1.5);
}

# Command line without extrema

    $cmd = qq[pcomb $pc0args $pc1args "$picture"];
    $cmd .= qq[ "$cpict"] if ($cpict);
    $cmd .= qq[ | pcompos -b 0 0 0 $scolpic 0 $sh0 +t .1 $slabinvpic 2 -1 ];
    $cmd .= qq[ -t .5 $slabpic 0 0 - $legwidth 0];

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
    $minval =~ s/(\.[0-9]{$scaledigits})[0-9]*/$1/;
    my $maxval = ($rmax * .27 + $gmax * .67 + $bmax * .06) * $mult;
    $maxval =~ s/(\.[0-9]{$scaledigits})[0-9]*/$1/;

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
