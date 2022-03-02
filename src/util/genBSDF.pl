#!/usr/bin/perl -w
# RCSid $Id: genBSDF.pl,v 2.89 2022/03/02 21:15:18 greg Exp $
#
# Compute BSDF based on geometry and material description
#
#	G. Ward
#
use strict;
my $windoz = ($^O eq "MSWin32" or $^O eq "MSWin64");
use File::Temp qw/ :mktemp  /;
sub userror {
	print STDERR "Usage: genBSDF [-n Nproc][-c Nsamp][-W][-t{3|4} Nlog2][-r \"ropts\"][-s \"x=string;y=string\"][-dim xmin xmax ymin ymax zmin zmax][{+|-}C][{+|-}a][{+|-}f][{+|-}b][{+|-}mgf][{+|-}geom units] [input ..]\n";
	exit 1;
}
my ($td,$radscn,$mgfscn,$octree,$fsender,$bsender,$receivers,$facedat,$behinddat,$rmtmp);
my ($tf,$rf,$tb,$rb,$tfx,$rfx,$tbx,$rbx,$tfz,$rfz,$tbz,$rbz,$cph);
my ($curphase, $recovery);
if ($#ARGV == 1 && "$ARGV[0]" =~ /^-rec/) {
	$td = $ARGV[1];
	open(MYAVH, "< $td/savedARGV.txt") or die "$td: invalid path\n";
	@ARGV = <MYAVH>;
	close MYAVH;
	chomp @ARGV;
	$recovery = 0;
	if (open(MYPH, "< $td/phase.txt")) {
		while (<MYPH>) {
			chomp($recovery = $_);
		}
		close MYPH;
	}
} elsif ($windoz) {
	my $tmploc = `echo \%TMP\%`;
	chomp $tmploc;
	$td = mkdtemp("$tmploc\\genBSDF.XXXXXX");
} else {
	$td = mkdtemp("/tmp/genBSDF.XXXXXX");
	chomp $td;
}
if ($windoz) {
	$radscn = "$td\\device.rad";
	$mgfscn = "$td\\device.mgf";
	$octree = "$td\\device.oct";
	$fsender = "$td\\fsender.rad";
	$bsender = "$td\\bsender.rad";
	$receivers = "$td\\receivers.rad";
	$facedat = "$td\\face.dat";
	$behinddat = "$td\\behind.dat";
	$tf = "$td\\tf.dat";
	$rf = "$td\\rf.dat";
	$tb = "$td\\tb.dat";
	$rb = "$td\\rb.dat";
	$tfx = "$td\\tfx.dat";
	$rfx = "$td\\rfx.dat";
	$tbx = "$td\\tbx.dat";
	$rbx = "$td\\rbx.dat";
	$tfz = "$td\\tfz.dat";
	$rfz = "$td\\rfz.dat";
	$tbz = "$td\\tbz.dat";
	$rbz = "$td\\rbz.dat";
	$cph = "$td\\phase.txt";
	$rmtmp = "rd /S /Q $td";
} else {
	$radscn = "$td/device.rad";
	$mgfscn = "$td/device.mgf";
	$octree = "$td/device.oct";
	$fsender = "$td/fsender.rad";
	$bsender = "$td/bsender.rad";
	$receivers = "$td/receivers.rad";
	$facedat = "$td/face.dat";
	$behinddat = "$td/behind.dat";
	$tf = "$td/tf.dat";
	$rf = "$td/rf.dat";
	$tb = "$td/tb.dat";
	$rb = "$td/rb.dat";
	$tfx = "$td/tfx.dat";
	$rfx = "$td/rfx.dat";
	$tbx = "$td/tbx.dat";
	$rbx = "$td/rbx.dat";
	$tfz = "$td/tfz.dat";
	$rfz = "$td/rfz.dat";
	$tbz = "$td/tbz.dat";
	$rbz = "$td/rbz.dat";
	$cph = "$td/phase.txt";
	$rmtmp = "rm -rf $td";
}
my @savedARGV = @ARGV;
my $rfluxmtx = "rfluxmtx -ab 5 -ad 700 -lw 3e-6 -w-";
my $wrapper = "wrapBSDF";
my $tensortree = 0;
my $ttlog2 = 4;
my $dorecip = 1;
my $nsamp = 2000;
my $mgfin = 0;
my $geout = 1;
my $nproc = 1;
my $docolor = 0;
my $doforw = 0;
my $doback = 1;
my $pctcull = 90;
my $gunit = "meter";
my $curspec = "Visible";
my @dim;
# Get options
while ($#ARGV >= 0) {
	if ("$ARGV[0]" =~ /^[-+]m/) {
		$mgfin = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" eq "-r") {
		$rfluxmtx .= " $ARGV[1]";
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^[-+]g/) {
		$geout = ("$ARGV[0]" =~ /^\+/);
		$gunit = $ARGV[1];
		if ($gunit !~ /^(?i)(meter|foot|inch|centimeter|millimeter)$/) {
			die "Illegal geometry unit '$gunit': must be meter, foot, inch, centimeter, or millimeter\n";
		}
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^[-+]C/) {
		$docolor = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]a/) {
		$dorecip = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]f/) {
		$doforw = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]b/) {
		$doback = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" eq "-t") {
		# Use value < 0 for rttree_reduce bypass
		$pctcull = $ARGV[1];
		if ($pctcull >= 100) {
			die "Illegal -t culling percentage, must be < 100\n";
		}
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-t[34]$/) {
		$tensortree = substr($ARGV[0], 2, 1);
		$ttlog2 = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-s") {
		$wrapper .= " -f \"$ARGV[1]\"";
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-W") {
		$wrapper .= " -W";
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
$wrapper .= $tensortree ? " -a t$tensortree" : " -a kf -c";
$wrapper .= " -u $gunit";
if ( !defined $recovery ) {
	# Issue warning for unhandled reciprocity case
	print STDERR "Warning: recommend both +forward and +backward with -t3\n" if
			($tensortree==3 && !($doforw && $doback));
	# Get scene description
	if ( $mgfin ) {
		system "mgf2rad @ARGV > $radscn";
		die "Could not load MGF input\n" if ( $? );
	} else {
		system "xform -e @ARGV > $radscn";
		die "Could not load Radiance input\n" if ( $? );
	}
}
if ( $#dim != 5 ) {
	@dim = split ' ', `getbbox -h -w $radscn`;
}
die "Device entirely inside room!\n" if ( $dim[4] >= 0 );
if ( $dim[5] > 1e-5 ) {
	print STDERR "Warning: Device extends into room\n";
} elsif ( $dim[5]*$dim[5] > .01*($dim[1]-$dim[0])*($dim[3]-$dim[2]) ) {
	print STDERR "Warning: Device far behind Z==0 plane\n";
}
# Assume Zmax==0 to derive thickness so pkgBSDF will work
$wrapper .= ' -f "t=' . (-$dim[4]) . ';w=' . ($dim[1] - $dim[0]) .
		';h=' . ($dim[3] - $dim[2]) . '"';
$wrapper .= " -g $mgfscn" if ( $geout );
# Calculate CIE (u',v') from Radiance RGB:
my $CIEuv =	'Xi=.5141*Ri+.3239*Gi+.1620*Bi;' .
		'Yi=.2651*Ri+.6701*Gi+.0648*Bi;' .
		'Zi=.0241*Ri+.1229*Gi+.8530*Bi;' .
		'den=Xi+15*Yi+3*Zi;' .
		'uprime=if(Yi,4*Xi/den,4/19);' .
		'vprime=if(Yi,9*Yi/den,9/19);' ;
my $FEPS = 1e-5;
my $ns = 2**$ttlog2;
my $nx = int(sqrt($nsamp*($dim[1]-$dim[0])/($dim[3]-$dim[2])) + 1);
my $ny = int($nsamp/$nx + 1);
$nsamp = $nx * $ny;
$rfluxmtx .= " -n $nproc -c $nsamp";
if ( !defined $recovery ) {
	open(MYAVH, "> $td/savedARGV.txt");
	foreach (@savedARGV) {
		print MYAVH "$_\n";
	}
	close MYAVH;
	# Generate octree
	system "oconv -w $radscn > $octree";
	die "Could not compile scene\n" if ( $? );
	# Add MGF description if requested
	if ( $geout ) {
		open(MGFSCN, "> $mgfscn");
		printf MGFSCN "xf -t %.6f %.6f 0\n", -($dim[0]+$dim[1])/2, -($dim[2]+$dim[3])/2;
		close MGFSCN;
		if ( $mgfin ) {
			die "+mgf requires input file with +geom\n" if ($#ARGV < 0);
			system qq{mgfilt "#,o,xf,c,cxy,cspec,cmix,m,sides,rd,td,rs,ts,ir,v,p,n,f,fh,sph,cyl,cone,prism,ring,torus" @ARGV >> $mgfscn};
		} else {
			system "rad2mgf $radscn >> $mgfscn";
		}
		open(MGFSCN, ">> $mgfscn");
		print MGFSCN "xf\n";
		close MGFSCN;
	}
	# Create receiver & sender surfaces (rectangular)
	open(RADSCN, "> $receivers");
	print RADSCN '#@rfluxmtx ' . ($tensortree ? "h=-sc$ns\n" : "h=-kf\n");
	print RADSCN '#@rfluxmtx ' . "u=-Y o=$facedat\n\n";
	print RADSCN "void glow receiver_face\n0\n0\n4 1 1 1 0\n\n";
	print RADSCN "receiver_face source f_receiver\n0\n0\n4 0 0 1 180\n\n";
	print RADSCN '#@rfluxmtx ' . ($tensortree ? "h=+sc$ns\n" : "h=+kf\n");
	print RADSCN '#@rfluxmtx ' . "u=-Y o=$behinddat\n\n";
	print RADSCN "void glow receiver_behind\n0\n0\n4 1 1 1 0\n\n";
	print RADSCN "receiver_behind source b_receiver\n0\n0\n4 0 0 -1 180\n";
	close RADSCN;
	# Prepare sender surfaces
	if ( $tensortree != 3 ) {	# Isotropic tensor tree is exception
		open (RADSCN, "> $fsender");
		print RADSCN '#@rfluxmtx u=-Y ' . ($tensortree ? "h=-sc$ns\n\n" : "h=-kf\n\n");
		print RADSCN "void polygon fwd_sender\n0\n0\n12\n";
		printf RADSCN "\t%e\t%e\t%e\n", $dim[0], $dim[2], $dim[4]-$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[0], $dim[3], $dim[4]-$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[1], $dim[3], $dim[4]-$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[1], $dim[2], $dim[4]-$FEPS;
		close RADSCN;
		open (RADSCN, "> $bsender");
		print RADSCN '#@rfluxmtx u=-Y ' . ($tensortree ? "h=+sc$ns\n\n" : "h=+kf\n\n");
		print RADSCN "void polygon bwd_sender\n0\n0\n12\n";
		printf RADSCN "\t%e\t%e\t%e\n", $dim[0], $dim[2], $dim[5]+$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[1], $dim[2], $dim[5]+$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[1], $dim[3], $dim[5]+$FEPS;
		printf RADSCN "\t%e\t%e\t%e\n", $dim[0], $dim[3], $dim[5]+$FEPS;
		close RADSCN;
	}
	print STDERR "Recover using: $0 -recover $td\n";
}
$curphase = 0;
# Create data segments (all the work happens here)
if ( $tensortree ) {
	do_tree_bsdf();
} else {
	do_matrix_bsdf();
}
# Output XML
# print STDERR "Running: $wrapper\n";
system "$wrapper -C \"Created by: genBSDF @savedARGV\"";
die "Could not wrap BSDF data\n" if ( $? );
# Clean up temporary files and exit
exec $rmtmp;

#============== End of main program segment ==============#

# Function to determine if next phase should be skipped or recovered
sub do_phase {
	$curphase++;
	if ( defined $recovery ) {
		if ( $recovery > $curphase ) { return 0; }
		if ( $recovery == $curphase ) { return -1; }
	}
	open(MYPH, ">> $td/phase.txt");
	print MYPH "$curphase\n";
	close MYPH;
	return 1;
}

# Check if we are in active phase (not skipping parts)
sub active_phase {
	if ( defined $recovery ) {
		if ( $recovery > $curphase ) { return 0; }
		if ( $recovery == $curphase ) { return -1; }
	}
	return 1;
}

# Function to run program and check output if in active phase
sub run_check {
	if ( !active_phase() ) { return; }
	my $cmd = shift;
	# print STDERR "Running: $cmd\n";
	system $cmd;
	die "Failure running: $cmd\n" if ( $? );
}

#++++++++++++++ Tensor tree BSDF generation ++++++++++++++#
sub do_tree_bsdf {

	# Run rfluxmtx processes to compute each side
	do_ttree_dir(0) if ( $doback );
	do_ttree_dir(1) if ( $doforw );

}	# end of sub do_tree_bsdf()

# Call rfluxmtx and process tensor tree BSDF for the given direction
sub do_ttree_dir {
	my $forw = shift;
	my $dop = do_phase();
	my $r = ($dop < 0) ? " -r" : "";
	my $cmd;
	if ( $tensortree == 3 ) {
		# Isotropic BSDF
		my $ns2 = $ns / 2;
		if ($windoz) {
			$cmd = "cnt $ns2 $ny $nx " .
				qq{| rcalc -e "r1=rand(.8681*recno-.673892)" } .
				qq{-e "r2=rand(-5.37138*recno+67.1737811)" } .
				qq{-e "r3=rand(+3.17603772*recno+83.766771)" } .
				qq{-e "r4=rand(-1.5839226*recno-59.82712)" } .
				qq{-e "odds(n):if(.5*n-floor(.5*n)-.25,-1,1)"} .
				qq{-e "Dx=1-(\$1+r1)/$ns2" } .
				qq{-e "Dy=min(1/$ns,sqrt(1-Dx*Dx))*odds(\$1)*r2" } .
				qq{-e "Dz=sqrt(1-Dx*Dx-Dy*Dy)" } .
				qq{-e "xp=(\$3+r2)*(($dim[1]-$dim[0])/$nx)+$dim[0]" } .
				qq{-e "yp=(\$2+r3)*(($dim[3]-$dim[2])/$ny)+$dim[2]" } .
				qq{-e "zp=$dim[5-$forw]" -e "myDz=Dz*($forw*2-1)" } .
				qq{-e "\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz" } .
				qq{-e "\$4=Dx;\$5=Dy;\$6=myDz" } .
				"| $rfluxmtx$r -fa -y $ns2 - $receivers -i $octree";
		} else {
			$cmd = "cnt $ns2 $ny $nx " .
				qq{| rcalc -e 'r1=rand(.8681*recno-.673892)' } .
				qq{-e 'r2=rand(-5.37138*recno+67.1737811)' } .
				qq{-e 'r3=rand(+3.17603772*recno+83.766771)' } .
				qq{-e 'r4=rand(-1.5839226*recno-59.82712)' } .
				qq{-e 'odds(n):if(.5*n-floor(.5*n)-.25,-1,1)' } .
				qq{-e 'Dx=1-(\$1+r1)/$ns2' } .
				qq{-e 'Dy=min(1/$ns,sqrt(1-Dx*Dx))*odds(\$1)*r2' } .
				qq{-e 'Dz=sqrt(1-Dx*Dx-Dy*Dy)' } .
				qq{-e 'xp=(\$3+r3)*(($dim[1]-$dim[0])/$nx)+$dim[0]' } .
				qq{-e 'yp=(\$2+r4)*(($dim[3]-$dim[2])/$ny)+$dim[2]' } .
				qq{-e 'zp=$dim[5-$forw]' -e 'myDz=Dz*($forw*2-1)' } .
				qq{-e '\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz' } .
				qq{-e '\$4=Dx;\$5=Dy;\$6=myDz' -of } .
				"| $rfluxmtx$r -h -ff -y $ns2 - $receivers -i $octree";
		}
	} else {
		# Anisotropic BSDF
		my $sender = ($bsender,$fsender)[$forw];
		if ($windoz) {
			$cmd = "$rfluxmtx$r -fa $sender $receivers -i $octree";
		} else {
			$cmd = "$rfluxmtx$r -h -ff $sender $receivers -i $octree";
		}
	}
	if ( $dop ) {
		# print STDERR "Running: $cmd\n";
		system $cmd;
		die "Failure running rfluxmtx" if ( $? );
	}
	ttree_out($forw);
}	# end of do_ttree_dir()

# Simplify and store tensor tree results
sub ttree_out {
	my $forw = shift;
	my ($refldat,$transdat);
	if ( $forw ) {
		$transdat = $facedat;
		$refldat = $behinddat;
	} else {
		$transdat = $behinddat;
		$refldat = $facedat;
	}
	# Only output one transmitted anisotropic distribution, preferring backwards
	if ( !$forw || !$doback || $tensortree==3 ) {
		my $ttyp = ("tb","tf")[$forw];
		ttree_comp($ttyp, "Visible", $transdat, ($tb,$tf)[$forw]);
		if ( $docolor ) {
			ttree_comp($ttyp, "CIE-u", $transdat, ($tbx,$tfx)[$forw]);
			ttree_comp($ttyp, "CIE-v", $transdat, ($tbz,$tfz)[$forw]);
		}
	}
	# Output reflection
	my $rtyp = ("rb","rf")[$forw];
	ttree_comp($rtyp, "Visible", $refldat, ($rb,$rf)[$forw]);
	if ( $docolor ) {
		ttree_comp($rtyp, "CIE-u", $refldat, ($rbx,$rfx)[$forw]);
		ttree_comp($rtyp, "CIE-v", $refldat, ($rbz,$rfz)[$forw]);
	}
}	# end of ttree_out()

# Call rttree_reduce on the given component
sub ttree_comp {
	my $typ = shift;
	my $spec = shift;
	my $src = shift;
	my $dest = shift;
	my $cmd;
	if ($windoz) {
		if ("$spec" eq "Visible") {
			$cmd = qq{rcalc -e "Omega:PI/($ns*$ns)" } .
				q{-e "Ri=$1;Gi=$2;Bi=$3" } .
				qq{-e "$CIEuv" } .
				q{-e "$1=Yi/Omega"};
		} elsif ("$spec" eq "CIE-u") {
			$cmd = q{rcalc -e "Ri=$1;Gi=$2;Bi=$3" } .
				qq{-e "$CIEuv" } .
				q{-e "$1=uprime"};
		} elsif ("$spec" eq "CIE-v") {
			$cmd = q{rcalc -e "Ri=$1;Gi=$2;Bi=$3" } .
				qq{-e "$CIEuv" } .
				q{-e "$1=vprime"};
		}
	} else {
		if ("$spec" eq "Visible") {
			$cmd = "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
				q{-e 'Ri=$1;Gi=$2;Bi=$3' } .
				"-e '$CIEuv' " .
				q{-e '$1=Yi/Omega'};
		} elsif ("$spec" eq "CIE-u") {
			$cmd = q{rcalc -if3 -e 'Ri=$1;Gi=$2;Bi=$3' } .
				"-e '$CIEuv' " .
				q{-e '$1=uprime'};
		} elsif ("$spec" eq "CIE-v") {
			$cmd = q{rcalc -if3 -e 'Ri=$1;Gi=$2;Bi=$3' } .
				"-e '$CIEuv' " .
				q{-e '$1=vprime'};
		}
	}
	if ($pctcull >= 0) {
		my $avg = ( $dorecip && ( $tensortree == 3 || "$typ" =~ /^r[fb]/ ) ) ? " -a" : "";
		my $pcull = ("$spec" eq "Visible") ? $pctcull :
						     (100 - (100-$pctcull)*.25) ;
		if ($windoz) {
			$cmd = "rcollate -ho -oc 1 $src | " .
					$cmd .
					" | rttree_reduce$avg -h -fa -t $pcull -r $tensortree -g $ttlog2";
		} else {
			$cmd .= " -of $src " .
					"| rttree_reduce$avg -h -ff -t $pcull -r $tensortree -g $ttlog2";
		}
		run_check "$cmd > $dest";
	} else {
		if ($windoz) {
			$cmd = "rcollate -ho -oc 1 $src | " . $cmd ;
		} else {
			$cmd .= " $src";
		}
		if ( active_phase() ) {
			open(DATOUT, "> $dest");
			print DATOUT "{\n";
			close DATOUT;
			# print STDERR "Running: $cmd\n";
			system "$cmd >> $dest";
			die "Failure running rcalc" if ( $? );
			open(DATOUT, ">> $dest");
			for (my $i = ($tensortree==3)*$ns*$ns*$ns/2; $i-- > 0; ) {
				print DATOUT "0\n";
			}
			print DATOUT "}\n";
			close DATOUT;
		}
	}
	if ( "$spec" ne "$curspec" ) {
		$wrapper .= " -s $spec";
		$curspec = $spec;
	}
	$wrapper .= " -$typ $dest";
}	# end of ttree_comp()

#------------- End of do_tree_bsdf() & subroutines -------------#

#+++++++++++++++ Klems matrix BSDF generation +++++++++++++++#
sub do_matrix_bsdf {

	# Run rfluxmtx processes to compute each side
	do_matrix_dir(0) if ( $doback );
	do_matrix_dir(1) if ( $doforw );

}	# end of sub do_matrix_bsdf()

# Call rfluxmtx and process tensor tree BSDF for the given direction
sub do_matrix_dir {
	my $forw = shift;
	my $dop = do_phase();
	my $r = ($dop < 0) ? " -r" : "";
	my $sender = ($bsender,$fsender)[$forw];
	my $cmd = "$rfluxmtx$r -fd $sender $receivers -i $octree";
	if ( $dop ) {
		# print STDERR "Running: $cmd\n";
		system $cmd;
		die "Failure running rfluxmtx" if ( $? );
	}
	matrix_out($forw);
}	# end of do_matrix_dir()

sub matrix_out {
	my $forw = shift;
	my ($refldat,$transdat);
	if ( $forw ) {
		$transdat = $facedat;
		$refldat = $behinddat;
	} else {
		$transdat = $behinddat;
		$refldat = $facedat;
	}
	# Output transmission
	my $ttyp = ("tb","tf")[$forw];
	matrix_comp($ttyp, "Visible", $transdat, ($tb,$tf)[$forw]);
	if ( $docolor ) {
		matrix_comp($ttyp, "CIE-X", $transdat, ($tbx,$tfx)[$forw]);
		matrix_comp($ttyp, "CIE-Z", $transdat, ($tbz,$tfz)[$forw]);
	}
	# Output reflection
	my $rtyp = ("rb","rf")[$forw];
	matrix_comp($rtyp, "Visible", $refldat, ($rb,$rf)[$forw]);
	if ( $docolor ) {
		matrix_comp($rtyp, "CIE-X", $refldat, ($rbx,$rfx)[$forw]);
		matrix_comp($rtyp, "CIE-Z", $refldat, ($rbz,$rfz)[$forw]);
	}
}	# end of matrix_out()

# Transpose matrix component data and save to file
sub matrix_comp {
	my $typ = shift;
	my $spec = shift;
	my $src = shift;
	my $dest = shift;
	my $cmd = "rmtxop -fa -t";
	if ("$spec" eq "Visible") {
		$cmd .= " -c 0.2651 0.6701 0.0648";
	} elsif ("$spec" eq "CIE-X") {
		$cmd .= " -c 0.5141 0.3239 0.1620";
	} elsif ("$spec" eq "CIE-Z") {
		$cmd .= " -c 0.0241 0.1229 0.8530";
	}
	$cmd .= " $src | getinfo -";
	run_check "$cmd > $dest";
	if ( "$spec" ne "$curspec" ) {
		$wrapper .= " -s $spec";
		$curspec = $spec;
	}
	$wrapper .= " -$typ $dest";
}	# end of matrix_comp()

#------------- End of do_matrix_bsdf() & subroutines --------------#
