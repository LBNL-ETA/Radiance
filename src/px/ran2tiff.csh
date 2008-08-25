#!/bin/csh -fe
# Convert Radiance animation frames to TIFF output
#
set histosiz=200
set pfwt=0.9
set outdir=""
set pcopts=()
set tfopts=()
if (! $#argv) set argv=(DUMMY)
# Process options for pcond and ra_tiff
while ("$argv[1]" =~ -*)
	switch ("$argv[1]")
	case -W:
		shift argv
		set pfwt=$argv[1]
		breaksw
	case -H:
		shift argv
		set histof=$argv[1]:q
		breaksw
	case -D:
		shift argv
		if (! -d $argv[1]:q ) then
			echo "Directory $argv[1] does not exist"
			exit 1
		endif
		set outdir=$argv[1]:q/
		breaksw
	case -h*:
	case -a*:
	case -v*:
	case -s*:
	case -c*:
	case -l*:
		set pcopts=($pcopts $argv[1])
		breaksw
	case -u:
	case -d:
	case -f:
		set pcopts=($pcopts $argv[1-2])
		shift argv
		breaksw
	case -p:
		shift argv
		set pcopts=($pcopts -p $argv[1-6])
		shift argv; shift argv; shift argv; shift argv; shift argv
		breaksw
	case -z:
	case -b:
	case -w:
		set tfopts=($tfopts $argv[1])
		breaksw
	case -g:
		shift argv
		set tfopts=($tfopts -g $argv[1])
		breaksw
	default:
		echo "$0: bad option: $argv[1]"
		exit 1
	endsw
	shift argv
end
if ($#argv < 2) then
	echo Usage: "$0 [-W prev_frame_wt][-H histo][-D dir][pcond opts][ra_tiff opts] frame1 frame2 .."
	exit 1
endif
# Get shrunken image luminances
set vald=`mktemp -d /tmp/val.XXXXXX`
foreach inp ($argv:q)
	set datf="$inp:t"
	set datf="$vald/$datf:r.dat"
	pfilt -1 -x 128 -y 128 -p 1 $inp:q \
		| pvalue -o -h -H -b -df \
		| rcalc -if1 -e 'L=$1*179;cond=L-1e-7;$1=log10(L)' \
		> $datf:q
end
# Get Min. and Max. log values
set Lmin=`cat $vald/*.dat | total -l | rcalc -e '$1=$1-.01'`
set Lmax=`cat $vald/*.dat | total -u | rcalc -e '$1=$1+.01'`
if ($?histof) then
	if (-r $histof) then
		# Fix min/max and translate histogram
		set Lmin=`sed -n '1p' $histof | rcalc -e 'min(a,b):if(a-b,b,a);$1=min($1,'"$Lmin)"`
		set Lmax=`sed -n '$p' $histof | rcalc -e 'max(a,b):if(a-b,a,b);$1=max($1,'"$Lmax)"`
		tabfunc -i hfunc < $histof > $vald/oldhist.cal
		cnt $histosiz \
			| rcalc -e "L10=$Lmin+($Lmax-$Lmin)/$histosiz"'*($1+.5)' \
				-f $vald/oldhist.cal -e '$1=L10;$2=hfunc(L10)' \
			> $vald/oldhisto.dat
	endif
endif
foreach inp ($argv:q)
	set datf="$inp:t"
	set datf="$vald/$datf:r.dat"
	set outp="$inp:t"
	set outp="$outdir$outp:r.tif"
	endif
	histo $Lmin $Lmax $histosiz < $datf > $vald/newhisto.dat
	if (-f $vald/oldhisto.dat) then
		rlam $vald/newhisto.dat $vald/oldhisto.dat \
			| rcalc -e '$1=$1;$2=$2+$4*'$pfwt \
			> $vald/histo.dat
	else
		mv $vald/{new,}histo.dat
	endif
	pcond $pcopts -I $inp:q < $vald/histo.dat \
		| ra_tiff $tfopts - $outp:q
	mv $vald/{,old}histo.dat
end
if ($?histof) then
	cp -f $vald/oldhisto.dat $histof
endif
rm -rf $vald
