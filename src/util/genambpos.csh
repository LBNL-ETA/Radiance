#!/bin/csh -f
# RCSid: $Id: genambpos.csh,v 3.1 2003/02/22 02:07:30 greg Exp $

# Mark ambient locations   J.Mardaljevic, edu@uk.ac.dmu

set ambdat=none
set radius=0.05

alias readvar 'echo -n Enter \!:1 "[$\!:1]: ";set ans="$<";if("$ans" != "")set \!:1="$ans"'

cat <<_EOF_

            CREATE SCENE DESCRIPTION MARKERS FOR
                 AMBIENT SAMPLING LOCATIONS

This script puts coloured markers (spheres) at points where hemispherical
sampling was initiated.  The locations and levels are taken from the binary
ambient file.  The levels are colour coded -

       Level 0 = red
       Level 1 = green
       Level 2 = blue
       Level 3 = yellow
       Level 4 = cyan
       Level 5 = magenta
       Level 6 = pink 

To view the locations, include the file "amb_pos.rad" when re-creating
the octree.

_EOF_

echo "Input ambient file"
readvar ambdat

echo "Input marker sphere radius"
readvar radius

getinfo $ambdat
lookamb -d -h $ambdat > af$$.tmp

set red=(1 0 0 1 0 1 1)
set grn=(0 1 0 1 1 0 0.5)
set blu=(0 0 1 0 1 1 0.5)

foreach j (1 2 3 4 5 6 7)

@ i = $j - 1

rcalc -e 'cond=eq($7,'"$i"');eq(a,b)=if(a-b+1e-7,b-a+1e-7,-1);$1=$1;$2=$2;$3=$3' af$$.tmp > amb_pos_lev$i

set empt=`wc -l < amb_pos_lev$i`

if ($empt !~ 0) then

echo "$empt ambient locations for level $i"

cat > sph_lev$i.fmt << _EOF_

mat_lev$i sphere sph_lev$i.\${recno}
0
0
4     \${ x1 } \${ y1 } \${ z1 }     $radius
_EOF_

cat >> amb_pos.rad << _EOF_
void plastic mat_lev$i
0
0
5   $red[$j]  $grn[$j]   $blu[$j]   0    0

!rcalc -e 'x1=\$1;y1=\$2;z1=\$3' -o sph_lev$i.fmt amb_pos_lev$i

_EOF_

else

rm amb_pos_lev$i

endif

end  

rm af$$.tmp
