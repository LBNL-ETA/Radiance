#!/bin/csh -f
# RCSid: $Id: compamb.csh,v 3.7 2008/08/25 04:50:32 greg Exp $
#
# Compute best ambient value for a scene and append to rad input file
#
while ( $#argv > 1 )
	switch ( $argv[1] )
	case -e:
		set doexpos
		breaksw
	case -c:
		set docolor
		breaksw
	default:
		goto userr
	endsw
	shift argv
end
userr:
if ( $#argv != 1 ) then
	echo Usage: $0 [-e][-c] rad_input_file
	exit 1
endif
onintr quit
set td=`mktemp -d /tmp/ca.XXXXXX`
set tf=$td/compamb
set oct=`rad -w -s -e -v 0 $argv[1] QUA=High AMB=$tf.amb OPT=$tf.opt | sed -n 's/^OCTREE= //p'`
rad -n -s -V $argv[1] \
	| rpict @$tf.opt -av 0 0 0 -aw 16 -dv- -S 1 -x 16 -y 16 -ps 1 $oct \
	| ra_rgbe - '\!pvalue -h -H -d' > $tf.dat
echo \# Rad input file modified by $0 `date` >> $argv[1]
if ( $?doexpos ) then
	(echo -n 'EXPOSURE= '; \
		total -u $tf.dat | rcalc -e '$1=2/(.265*$1+.670*$2+.065*$3)') \
		>> $argv[1]
endif
lookamb -h -d $tf.amb | rcalc -e '$1=$10;$2=$11;$3=$12' >> $tf.dat
set lavg=`rcalc -e '$1=lum;lum=.265*$1+.670*$2+.065*$3;cond=lum-1e-5' $tf.dat | total -m -p`
if ( $?docolor ) then
	set cavg=(`total -m $tf.dat`)
	set av=(`rcalc -n -e "r=$cavg[1];g=$cavg[2];b=$cavg[3];sf=$lavg/(.265*r+.670*g+.065*b)" -e '$1=sf*r;$2=sf*g;$3=sf*b'`)
else
	set av=($lavg $lavg $lavg)
endif
echo "render= -av $av" >> $argv[1]
quit:
exec rm -r $td
