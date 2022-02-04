#!/usr/bin/perl -w
# RCSid $Id: phisto.pl,v 2.2 2022/02/04 20:11:49 greg Exp $
#
# Compute foveal histogram for picture set
#
#	G. Ward
#
use strict;
my $windoz = ($^O eq "MSWin32" or $^O eq "MSWin64");
use File::Temp qw/ :mktemp  /;
my $minop = 'L=$1*179;$1=if(L-1e-7,log10(L)-.01,-7)';
my $maxop = '$1=log10($1*179)+.01';
my $tf;
if ($windoz) {
	my $tmploc = `echo \%TMP\%`;
	chomp $tmploc;
	$tf = mktemp("$tmploc\\phdXXXXXX");
	chomp $tf;
	if ($#ARGV < 0) {
		system "pfilt -1 -x 128 -y 128 -p 1 " .
			"| pvalue -O -h -H -d -b > $tf";
	} else {
		foreach (@ARGV) {
			system "pfilt -1 -x 128 -y 128 -p 1 $_" .
				"| pvalue -O -h -H -d -b >> $tf";
			die "Bad picture '$_'\n" if ( $? );
		}
	}
	my $Lmin=`total -l $tf | rcalc -e "$minop"`;
	my $Lmax=`total -u $tf | rcalc -e "$maxop"`;
	chomp $Lmin;
	chomp $Lmax;
	system q[rcalc -e "L=$1*179;cond=L-1e-7;$1=log10(L)" ] . $tf .
		" | histo $Lmin $Lmax 777";
} else {
	$tf = mktemp("/tmp/phdXXXXXX");
	chomp $tf;
	if ($#ARGV < 0) {
		system "pfilt -1 -x 128 -y 128 -p 1 " .
			"| pvalue -O -h -H -df -b > $tf";
	} else {
		foreach (@ARGV) {
			system "pfilt -1 -x 128 -y 128 -p 1 '$_'" .
				"| pvalue -O -h -H -df -b >> $tf";
			die "Bad picture '$_'\n" if ( $? );
		}
	}
	my $Lmin=`total -if -l $tf | rcalc -e '$minop'`;
	my $Lmax=`total -if -u $tf | rcalc -e '$maxop'`;
	chomp $Lmin;
	chomp $Lmax;
	system q[rcalc -if -e 'L=$1*179;cond=L-1e-7;$1=log10(L)' ] . $tf .
		" | histo $Lmin $Lmax 777";
}
unlink $tf;
