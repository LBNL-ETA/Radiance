#!/usr/bin/perl -w
# RCSid $Id: genBSDF.pl,v 2.51 2015/03/12 13:24:41 greg Exp $
#
# Compute BSDF based on geometry and material description
#
#	G. Ward
#
use strict;
my $windoz = ($^O eq "MSWin32" or $^O eq "MSWin64");
use File::Temp qw/ :mktemp  /;
sub userror {
	print STDERR "Usage: genBSDF [-n Nproc][-c Nsamp][-W][-t{3|4} Nlog2][-r \"ropts\"][-f \"x=string;y=string\"][-dim xmin xmax ymin ymax zmin zmax][{+|-}C][{+|-}f][{+|-}b][{+|-}mgf][{+|-}geom units] [input ..]\n";
	exit 1;
}
my ($td,$radscn,$mgfscn,$octree,$fsender,$bsender,$receivers,$facedat,$behinddat,$rmtmp);
my ($tf,$rf,$tb,$rb,$tfx,$rfx,$tbx,$rbx,$tfz,$rfz,$tbz,$rbz);
if ($windoz) {
	my $tmploc = `echo \%TMP\%`;
	chomp($tmploc);
	$td = mkdtemp("$tmploc\\genBSDF.XXXXXX");
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
	chomp $td;
	$rmtmp = "rmdir /S /Q $td";
} else {
	$td = mkdtemp("/tmp/genBSDF.XXXXXX");
	chomp $td;
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
	$rmtmp = "rm -rf $td";
}
my @savedARGV = @ARGV;
my $rfluxmtx = "rfluxmtx -v -ab 5 -ad 700 -lw 3e-6";
my $wrapper = "wrapBSDF";
my $tensortree = 0;
my $ttlog2 = 4;
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
	} elsif ("$ARGV[0]" =~ /^[-+]f/) {
		$doforw = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" =~ /^[-+]b/) {
		$doback = ("$ARGV[0]" =~ /^\+/);
	} elsif ("$ARGV[0]" eq "-t") {
		# Use value < 0 for rttree_reduce bypass
		$pctcull = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" =~ /^-t[34]$/) {
		$tensortree = substr($ARGV[0], 2, 1);
		$ttlog2 = $ARGV[1];
		shift @ARGV;
	} elsif ("$ARGV[0]" eq "-f") {
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
# Issue warning for unhandled reciprocity case
print STDERR "Warning: recommend both +forward and +backward with -t3\n" if
		($tensortree==3 && !($doforw && $doback));
$wrapper .= $tensortree ? " -a t$tensortree" : " -a kf -c";
$wrapper .= " -u $gunit";
# Get scene description and dimensions
if ( $mgfin ) {
	system "mgf2rad @ARGV > $radscn";
	die "Could not load MGF input\n" if ( $? );
} else {
	system "xform -e @ARGV > $radscn";
	die "Could not load Radiance input\n" if ( $? );
}
if ($#dim != 5) {
	@dim = split ' ', `getbbox -h $radscn`;
}
print STDERR "Warning: Device extends into room\n" if ($dim[5] > 1e-5);
$wrapper .= ' -f "t=' . ($dim[5] - $dim[4]) . ';w=' . ($dim[1] - $dim[0]) .
		';h=' . ($dim[3] - $dim[2]) . '"';
# Generate octree
system "oconv -w $radscn > $octree";
die "Could not compile scene\n" if ( $? );
# Add MGF description if requested
if ( $geout ) {
	open(MGFSCN, "> $mgfscn");
	printf MGFSCN "xf -t %.6f %.6f 0\n", -($dim[0]+$dim[1])/2, -($dim[2]+$dim[3])/2;
	close MGFSCN;
	if ( $mgfin ) {
		system qq{mgfilt "#,o,xf,c,cxy,cspec,cmix,m,sides,rd,td,rs,ts,ir,v,p,n,f,fh,sph,cyl,cone,prism,ring,torus" @ARGV >> $mgfscn};
	} else {
		system "rad2mgf $radscn >> $mgfscn";
	}
	open(MGFSCN, ">> $mgfscn");
	print MGFSCN "xf\n";
	close MGFSCN;
	$wrapper .= " -g $mgfscn";
}
# Create receiver & sender surfaces (rectangular)
my $nx = int(sqrt($nsamp*($dim[1]-$dim[0])/($dim[3]-$dim[2])) + 1);
my $ny = int($nsamp/$nx + 1);
$nsamp = $nx * $ny;
my $ns = 2**$ttlog2;
open(RADSCN, "> $receivers");
print RADSCN '#@rfluxmtx ' . ($tensortree ? "h=sc$ns\n" : "h=kf\n");
print RADSCN '#@rfluxmtx ' . "u=Y o=$facedat\n\n";
print RADSCN "void glow receiver_face\n0\n0\n4 1 1 1 0\n\n";
print RADSCN "receiver_face source f_receiver\n0\n0\n4 0 0 1 180\n\n";
print RADSCN '#@rfluxmtx ' . "u=Y o=$behinddat\n\n";
print RADSCN "void glow receiver_behind\n0\n0\n4 1 1 1 0\n\n";
print RADSCN "receiver_behind source b_receiver\n0\n0\n4 0 0 -1 180\n";
close RADSCN;
# Finish rfluxmtx command and prepare sender surfaces
$rfluxmtx .= " -n $nproc -c $nsamp";
if ( $tensortree != 3 ) {	# Isotropic tensor tree is exception
	open (RADSCN, "> $fsender");
	print RADSCN '#@rfluxmtx u=Y ' . ($tensortree ? "h=sc$ns\n\n" : "h=kf\n\n");
	print RADSCN "void polygon fwd_sender\n0\n0\n12\n";
	printf RADSCN "\t%f\t%f\t%f\n", $dim[0], $dim[2], $dim[4];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[0], $dim[3], $dim[4];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[1], $dim[3], $dim[4];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[1], $dim[2], $dim[4];
	close RADSCN;
	open (RADSCN, "> $bsender");
	print RADSCN '#@rfluxmtx u=Y ' . ($tensortree ? "h=sc$ns\n\n" : "h=kf\n\n");
	print RADSCN "void polygon fwd_sender\n0\n0\n12\n";
	printf RADSCN "\t%f\t%f\t%f\n", $dim[0], $dim[2], $dim[5];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[1], $dim[2], $dim[5];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[1], $dim[3], $dim[5];
	printf RADSCN "\t%f\t%f\t%f\n", $dim[0], $dim[3], $dim[5];
	close RADSCN;
}
# Create data segments (all the work happens here)
if ( $tensortree ) {
	do_tree_bsdf();
} else {
	do_matrix_bsdf();
}
# Output XML
my $old_fh = select(STDOUT);
$| = 1;
select($old_fh);
print "<!-- Created by: genBSDF @savedARGV -->\n";
# print STDERR "Running: $wrapper\n";
system $wrapper;
die "Could not wrap BSDF data\n" if ( $? );
# Clean up temporary files and exit
system $rmtmp;

#-------------- End of main program segment --------------#

#++++++++++++++ Tensor tree BSDF generation ++++++++++++++#
sub do_tree_bsdf {

	# Run rfluxmtx processes to compute each side
	do_ttree_dir(0) if ( $doback );
	do_ttree_dir(1) if ( $doforw );

}	# end of sub do_tree_bsdf()

# Call rfluxmtx and process tensor tree BSDF for the given direction
sub do_ttree_dir {
	my $forw = shift;
	my $cmd;
	if ( $tensortree == 3 ) {
		# Isotropic BSDF
		my $ns2 = $ns / 2;
		if ($windoz) {
			$cmd = "cnt $ns2 $ny $nx " .
				qq{| rcalc -e "r1=rand(.8681*recno-.673892)" } .
				qq{-e "r2=rand(-5.37138*recno+67.1737811)" } .
				qq{-e "r3=rand(+3.17603772*recno+83.766771)" } .
				qq{-e "Dx=1-2*(\$1+r1)/$ns;Dy:0;Dz=sqrt(1-Dx*Dx)" } .
				qq{-e "xp=(\$3+r2)*(($dim[1]-$dim[0])/$nx)+$dim[0]" } .
				qq{-e "yp=(\$2+r3)*(($dim[3]-$dim[2])/$ny)+$dim[2]" } .
				qq{-e "zp=$dim[5-$forw]" -e "myDz=Dz*($forw*2-1)" } .
				qq{-e "\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz" } .
				qq{-e "\$4=Dx;\$5=Dy;\$6=myDz" } .
				"| $rfluxmtx -fa -y $ns2 - $receivers -i $octree";
		} else {
			$cmd = "cnt $ns2 $ny $nx " .
				qq{| rcalc -e "r1=rand(.8681*recno-.673892)" } .
				qq{-e "r2=rand(-5.37138*recno+67.1737811)" } .
				qq{-e "r3=rand(+3.17603772*recno+83.766771)" } .
				qq{-e "Dx=1-2*(\$1+r1)/$ns;Dy:0;Dz=sqrt(1-Dx*Dx)" } .
				qq{-e "xp=(\$3+r2)*(($dim[1]-$dim[0])/$nx)+$dim[0]" } .
				qq{-e "yp=(\$2+r3)*(($dim[3]-$dim[2])/$ny)+$dim[2]" } .
				qq{-e "zp=$dim[5-$forw]" -e "myDz=Dz*($forw*2-1)" } .
				qq{-e '\$1=xp-Dx;\$2=yp-Dy;\$3=zp-myDz' } .
				qq{-e '\$4=Dx;\$5=Dy;\$6=myDz' -of } .
				"| $rfluxmtx -h -ff -y $ns2 - $receivers -i $octree";
		}
	} else {
		# Anisotropic BSDF
		my $sender = ($bsender,$fsender)[$forw];
		if ($windoz) {
			$cmd = "$rfluxmtx -fa $sender $receivers -i $octree";
		} else {
			$cmd = "$rfluxmtx -h -ff $sender $receivers -i $octree";
		}
	}
	# print STDERR "Starting: $cmd\n";
	system $cmd;
	die "Failure running rfluxmtx" if ( $? );
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
			ttree_comp($ttyp, "CIE-X", $transdat, ($tbx,$tfx)[$forw]);
			ttree_comp($ttyp, "CIE-Z", $transdat, ($tbz,$tfz)[$forw]);
		}
	}
	# Output reflection
	my $rtyp = ("rb","rf")[$forw];
	ttree_comp($rtyp, "Visible", $refldat, ($rb,$rf)[$forw]);
	if ( $docolor ) {
		ttree_comp($rtyp, "CIE-X", $refldat, ($rbx,$rfx)[$forw]);
		ttree_comp($rtyp, "CIE-Z", $refldat, ($rbz,$rfz)[$forw]);
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
				q{-e "$1=(0.265*$1+0.670*$2+0.065*$3)/Omega"};
		} elsif ("$spec" eq "CIE-X") {
			$cmd = qq{rcalc -e "Omega:PI/($ns*$ns)" } .
				q{-e "$1=(0.514*$1+0.324*$2+0.162*$3)/Omega"};
		} elsif ("$spec" eq "CIE-Z") {
			$cmd = qq{rcalc -e "Omega:PI/($ns*$ns)" } .
				q{-e "$1=(0.024*$1+0.123*$2+0.853*$3)/Omega"};
		}
	} else {
		if ("$spec" eq "Visible") {
			$cmd = "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
				q{-e '$1=(0.265*$1+0.670*$2+0.065*$3)/Omega'};
		} elsif ("$spec" eq "CIE-X") {
			$cmd = "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
				q{-e '$1=(0.514*$1+0.324*$2+0.162*$3)/Omega'};
		} elsif ("$spec" eq "CIE-Z") {
			$cmd = "rcalc -if3 -e 'Omega:PI/($ns*$ns)' " .
				q{-e '$1=(0.024*$1+0.123*$2+0.853*$3)/Omega'};
		}
	}
	if ($pctcull >= 0) {
		my $avg = ( "$typ" =~ /^r[fb]/ ) ? " -a" : "";
		if ($windoz) {
			$cmd = "rcollate -ho -oc 1 $src | " .
					$cmd .
					" | rttree_reduce$avg -h -fa -t $pctcull -r $tensortree -g $ttlog2";
		} else {
			$cmd .= " -of $src " .
					"| rttree_reduce$avg -h -ff -t $pctcull -r $tensortree -g $ttlog2";
		}
		# print STDERR "Running: $cmd\n";
		system "$cmd > $dest";
		die "Failure running rttree_reduce" if ( $? );
	} else {
		if ($windoz) {
			$cmd = "rcollate -ho -oc 1 $src | " . $cmd ;
		} else {
			$cmd .= " $src";
		}
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
	my $cmd;
	my $sender = ($bsender,$fsender)[$forw];
	$cmd = "$rfluxmtx -fd $sender $receivers -i $octree";
	# print STDERR "Starting: $cmd\n";
	system $cmd;
	die "Failure running rfluxmtx" if ( $? );
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
		$cmd .= " -c 0.265 0.670 0.065";
	} elsif ("$spec" eq "CIE-X") {
		$cmd .= " -c 0.514 0.324 0.162";
	} elsif ("$spec" eq "CIE-Z") {
		$cmd .= " -c 0.024 0.123 0.853";
	}
	$cmd .= " $src | rcollate -ho";
	# print STDERR "Running: $cmd\n";
	system "$cmd > $dest";
	die "Failure running rttree_reduce" if ( $? );
	if ( "$spec" ne "$curspec" ) {
		$wrapper .= " -s $spec";
		$curspec = $spec;
	}
	$wrapper .= " -$typ $dest";
}	# end of matrix_comp()

#------------- End of do_matrix_bsdf() & subroutines --------------#
