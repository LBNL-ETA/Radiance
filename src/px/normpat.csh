#!/bin/csh -f
# RCSid: $Id: normpat.csh,v 2.5 2004/01/29 22:18:26 greg Exp $
#
# Normalize a pattern for tiling (-b option blends edges) by removing
# lowest frequencies from image (-f option) and reducing to
# standard size (-r option)
#
set ha=$0
set ha=$ha:t
set pf="pfilt -e 2"
while ($#argv > 0)
	switch ($argv[1])
	case -r:
		shift argv
		set pf="$pf -x $argv[1] -y $argv[1] -p 1 -r .67"
		breaksw
	case -f:
		set ha="$ha -f"
		set dofsub
		breaksw
	case -b:
		set ha="$ha -b"
		set blend
		breaksw
	case -v:
		set verb
		breaksw
	case -*:
		echo bad option $argv[1]
		exit 1
	default:
		goto dofiles
	endsw
	shift argv
end
dofiles:
onintr quit
set td=/usr/tmp/np$$
mkdir $td
goto skipthis
cat > $td/coef.fmt << '_EOF_'
   rm:${  $25   };    gm:${  $26   };    bm:${  $27   };
  rcx:${   $1   };   gcx:${   $9   };   bcx:${  $17   };
  rcy:${   $2   };   gcy:${  $10   };   bcy:${  $18   };
  rsx:${   $3   };   gsx:${  $11   };   bsx:${  $19   };
  rsy:${   $4   };   gsy:${  $12   };   bsy:${  $20   };
rcxcy:${   $5   }; gcxcy:${  $13   }; bcxcy:${  $21   };
rcxsy:${   $6   }; gcxsy:${  $14   }; bcxsy:${  $22   };
rsxcy:${   $7   }; gsxcy:${  $15   }; bsxcy:${  $23   };
rsxsy:${   $8   }; gsxsy:${  $16   }; bsxsy:${  $24   };
'_EOF_'
cat > $td/coef.cal << '_EOF_'
$1=$3*2*cx; $2=$3*2*cy; $3=$3*2*sx; $4=$3*2*sy;
$5=$3*4*cx*cy; $6=$3*4*cx*sy; $7=$3*4*sx*cy; $8=$3*4*sx*sy;
$9=$4*2*cx; $10=$4*2*cy; $11=$4*2*sx; $12=$4*2*sy;
$13=$4*4*cx*cy; $14=$4*4*cx*sy; $15=$4*4*sx*cy; $16=$4*4*sx*sy;
$17=$5*2*cx; $18=$5*2*cy; $19=$5*2*sx; $20=$5*2*sy;
$21=$5*4*cx*cy; $22=$5*4*cx*sy; $23=$5*4*sx*cy; $24=$5*4*sx*sy;
$25=$3; $26=$4; $27=$5;
cx=cos(wx); cy=cos(wy);
sx=sin(wx); sy=sin(wy);
wx=2*PI/xres*($1+.5); wy=2*PI/yres*($2+.5);
'_EOF_'
cat > $td/fsub.cal << '_EOF_'
PI:3.14159265358979323846;
ro=ri(1)*rm/(rm+rcx*cx+rcy*cy+rsx*sx+rsy*sy
		+rcxcy*cx*cy+rcxsy*cx*sy+rsxcy*sx*cy+rsxsy*sx*sy);
go=gi(1)*gm/(gm+gcx*cx+gcy*cy+gsx*sx+gsy*sy
		+gcxcy*cx*cy+gcxsy*cx*sy+gsxcy*sx*cy+gsxsy*sx*sy);
bo=bi(1)*bm/(bm+bcx*cx+bcy*cy+bsx*sx+bsy*sy
		+bcxcy*cx*cy+bcxsy*cx*sy+bsxcy*sx*cy+bsxsy*sx*sy);
cx=cos(wx); cy=cos(wy);
sx=sin(wx); sy=sin(wy);
wx=2*PI/xres*(x+.5); wy=2*PI/yres*(y+.5);
'_EOF_'
skipthis:
foreach f ($*)
	if ( $?verb ) then
		echo $f\:
		echo adjusting exposure/size...
	endif
	$pf $f > $td/pf
	getinfo < $td/pf > $f
	ed - $f << _EOF_
i
$ha
.
w
q
_EOF_
	set resolu=`getinfo -d < $td/pf | sed 's/-Y \([1-9][0-9]*\) +X \([1-9][0-9]*\)/\2 \1/'`
	if ( ! $?dofsub ) then
		mv $td/pf $td/hf
		goto donefsub
	endif
	# if ( $?verb ) then
	#	echo computing Fourier coefficients...
	# endif
	# pfilt -1 -x 32 -y 32 $td/pf | pvalue -h \
	#	| rcalc -e 'xres:32;yres:32' -f $td/coef.cal \
	#	| total -m | rcalc -o $td/coef.fmt \
	#	> $td/coef
	if ( $?verb ) then
		# cat $td/coef
		echo removing low frequencies...
	endif
	pgblur -r `ev "sqrt($resolu[1]*$resolu[2])/8"` $td/pf > $td/lf
	pcomb -e 's=1/li(2);ro=s*ri(1);go=s*gi(1);bo=s*bi(1)' \
			$td/pf $td/lf > $td/hf
	# pcomb -f $td/coef -f $td/fsub.cal $td/pf > $td/hf
	donefsub:
	if ( $?blend ) then
		if ( $?verb ) then
			echo blending edges...
		endif
		@ mar= $resolu[1] - 3
		@ les= $resolu[1] - 1
		pcompos -x 1 $td/hf 0 0 | pfilt -1 -x 3 > $td/left
		pcompos $td/hf -$les 0 | pfilt -1 -x 3 > $td/right
		pcomb -e 'ro=f(ri);go=f(gi);bo=f(bi)' \
			-e 'f(p)=(3-x)/7*p(1)+(4+x)/7*p(2)' \
			$td/right $td/left > $td/left.patch
		pcomb -e 'ro=f(ri);go=f(gi);bo=f(bi)' \
			-e 'f(p)=(1+x)/7*p(1)+(6-x)/7*p(2)' \
			$td/left $td/right > $td/right.patch
		pcompos $td/hf 0 0 $td/left.patch 0 0 $td/right.patch $mar 0 \
			> $td/hflr
		@ mar= $resolu[2] - 3
		@ les= $resolu[2] - 1
		pcompos -y 1 $td/hflr 0 0 | pfilt -1 -y 3 > $td/bottom
		pcompos $td/hflr 0 -$les | pfilt -1 -y 3 > $td/top
		pcomb -e 'ro=f(ri);go=f(gi);bo=f(bi)' \
			-e 'f(p)=(3-y)/7*p(1)+(4+y)/7*p(2)' \
			$td/top $td/bottom > $td/bottom.patch
		pcomb -e 'ro=f(ri);go=f(gi);bo=f(bi)' \
			-e 'f(p)=(1+y)/7*p(1)+(6-y)/7*p(2)' \
			$td/bottom $td/top > $td/top.patch
		pcompos $td/hflr 0 0 $td/bottom.patch 0 0 $td/top.patch 0 $mar \
			| getinfo - >> $f
	else
		getinfo - < $td/hf >> $f
	endif
	if ( $?verb ) then
		echo $f done.
	endif
end
quit:
rm -rf $td
