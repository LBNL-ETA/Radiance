#!/bin/csh -f
# RCSid: $Id: glaze.csh,v 2.4 2004/08/25 01:40:12 greg Exp $
#
# Complex glazing model (goes with glaze1.cal and glaze2.cal)
#
# Oct. 2002 Greg Ward
# Aug. 2004 GW (added -f option to read glazings from file)
# Funding for this development generously provided by Visarc, Inc.
# (http://www.visarc.com)
#

#################################################################
#
# Supported surface types:
#
set sn_arr=("clear glass" "VE1-2M low-E coating" "PVB laminated" "V-175 white frit" "V-933 warm gray frit")
# Glass-side hemispherical reflectances for each surface type:
set rg_r_arr=(0.074 0.065 .11 0.33 0.15)
set rg_g_arr=(0.077 0.058 .11 0.33 0.15)
set rg_b_arr=(0.079 0.067 .11 0.33 0.15)
# Coating-side hemispherical reflectance for each surface type:
set rc_r_arr=(0.074 0.042 .11 0.59 0.21)
set rc_g_arr=(0.077 0.049 .11 0.59 0.21)
set rc_b_arr=(0.079 0.043 .11 0.59 0.21)
# Hemispherical (normal) transmittance for each surface type:
set tn_r_arr=(0.862 0.756 0.63 0.21 0.09)
set tn_g_arr=(0.890 0.808 0.63 0.21 0.09)
set tn_b_arr=(0.886 0.744 0.63 0.21 0.09)
# Boolean whether coatings can have partial coverage:
set part_arr=(0 0 0 1 1)

while ($#argv > 0)
	set header="Surface	Tr	Tg	Tb	Rcr	Rcg	Rcb	Rgr	Rgg	Rgb	Part"
	if ($#argv < 2 || "$argv[1]" != '-f') then
		echo "Usage: $0 [-f glazing.dat ..]"
		exit 1
	endif
	shift argv
	set gf="$argv[1]"
	shift argv
	if ("`sed -n 1p $gf:q`" != "$header") then
		echo "Bad header in $gf -- Expected: $header"
		exit 1
	endif
	echo "Adding glazing types from file $gf :"
	set nl=`wc -l < $gf:q`
	@ i=2
	while ($i <= $nl)
		set ln=(`sed -n ${i}p $gf:q`)
		if ($#ln != 11) then
			echo "Expected 11 words in line: $ln"
			exit 1
		endif
		echo $ln[1]
		set sn_arr=($sn_arr:q $ln[1])
		set tn_r_arr=($tn_r_arr $ln[2])
		set tn_g_arr=($tn_g_arr $ln[3])
		set tn_b_arr=($tn_b_arr $ln[4])
		set rc_r_arr=($rc_r_arr $ln[5])
		set rc_g_arr=($rc_g_arr $ln[6])
		set rc_b_arr=($rc_b_arr $ln[7])
		set rg_r_arr=($rc_r_arr $ln[8])
		set rg_g_arr=($rc_g_arr $ln[9])
		set rg_b_arr=($rc_b_arr $ln[10])
		set part_arr=($part_arr $ln[11])
		@ i++
	end
end

#################################################################
#
# Get user input
#
echo -n "Enter the number of panes in the system: "
set np="$<"
if ($np != 1 && $np != 2) then
	echo "Must be 1 or 2 pane system"
	exit 1
endif
echo ""
echo "Window normal faces interior"
echo ""
if ($np == 1) then
echo "  | |"
echo "  | |"
echo "  | |"
echo "  | |-->"
echo "  | |"
echo "  | |"
echo "  | |"
echo " s1 s2"
else
echo "  | |    | |"
echo "  | |    | |"
echo "  | |    | |"
echo "  | |    | |-->"
echo "  | |    | |"
echo "  | |    | |"
echo "  | |    | |"
echo " s1 s2  s3 s4"
endif
echo ""
echo "Supported surface types are:"
set i=1
while ($i <= $#sn_arr)
	echo "  " $i - $sn_arr[$i]
	@ i++
end
echo ""
echo -n "What is the type of s1? "
set s1t="$<"
echo -n "What is the type of s2? "
set s2t="$<"
if ($s1t != 1 && $s2t != 1) then
	echo "One surface of each pane must be $sn_arr[1]"
	exit 1
endif
if ($part_arr[$s1t]) then
	echo -n "Enter fraction coverage for s1 (0-1): "
	set s1c="$<"
endif
if ($part_arr[$s2t]) then
	echo -n "Enter fraction coverage for s2 (0-1): "
	set s2c="$<"
endif
if ($np == 2) then
echo -n "What is the type of s3? "
set s3t="$<"
echo -n "What is the type of s4? "
set s4t="$<"
if ($s3t != 1 && $s4t != 1) then
	echo "One surface of each pane must be $sn_arr[1]"
	exit 1
endif
if ($part_arr[$s3t]) then
	echo -n "Enter fraction coverage for s3 (0-1): "
	set s3c="$<"
endif
if ($part_arr[$s4t]) then
	echo -n "Enter fraction coverage for s4 (0-1): "
	set s4c="$<"
endif
endif

#################################################################
#
# Begin material comments
#
echo ""
echo "############################################"
echo "# Glazing produced by Radiance glaze script"
echo "# `date`"
echo "# Material surface normal points to interior"
echo "# Number of panes in system: $np"

if ($np == 2) goto glaze2
#################################################################
#
# Compute single glazing
#
set sc=1
echo "# Exterior surface s1 type: $sn_arr[$s1t]"
if ($?s1c) then
	echo "# s1 coating coverage: $s1c"
	set sc=$s1c
endif
echo "# Interior surface s2 type: $sn_arr[$s2t]"
if ($?s2c) then
	echo "# s2 coating coverage: $s2c"
	set sc=$s2c
endif
if ($s1t != 1) then
	set ct=$s1t
	echo -n "# Exterior normal hemispherical reflectance: "
	ev ".265*($sc*$rc_r_arr[$ct]+(1-$sc)*$rc_r_arr[1])+.670*($sc*$rc_g_arr[$ct]+(1-$sc)*$rc_g_arr[1])+.065*($sc*$rc_b_arr[$ct]+(1-$sc)*$rc_b_arr[1])"
	echo -n "# Interior normal hemispherical reflectance: "
	ev ".265*($sc*$rg_r_arr[$ct]+(1-$sc)*$rg_r_arr[1])+.670*($sc*$rg_g_arr[$ct]+(1-$sc)*$rg_g_arr[1])+.065*($sc*$rg_b_arr[$ct]+(1-$sc)*$rg_b_arr[1])"
else
	set ct=$s2t
	echo -n "# Exterior normal hemispherical reflectance: "
	ev ".265*($sc*$rg_r_arr[$ct]+(1-$sc)*$rg_r_arr[1])+.670*($sc*$rg_g_arr[$ct]+(1-$sc)*$rg_g_arr[1])+.065*($sc*$rg_b_arr[$ct]+(1-$sc)*$rg_b_arr[1])"
	echo -n "# Interior normal hemispherical reflectance: "
	ev ".265*($sc*$rc_r_arr[$ct]+(1-$sc)*$rc_r_arr[1])+.670*($sc*$rc_g_arr[$ct]+(1-$sc)*$rc_g_arr[1])+.065*($sc*$rc_b_arr[$ct]+(1-$sc)*$rc_b_arr[1])"
endif
echo -n "# Normal hemispherical transmittance: "
ev ".265*($sc*$tn_r_arr[$ct]+(1-$sc)*$tn_r_arr[1])+.670*($sc*$tn_g_arr[$ct]+(1-$sc)*$tn_g_arr[1])+.065*($sc*$tn_b_arr[$ct]+(1-$sc)*$tn_b_arr[1])"
echo "#"
echo "void BRTDfunc glaze1_unnamed"
if ($part_arr[$s1t] || $part_arr[$s2t]) then
### Frit glazing
echo "10"
echo "	sr_frit_r sr_frit_g sr_frit_b"
echo "	st_frit_r st_frit_g st_frit_b"
echo "	0 0 0"
echo "	glaze1.cal"
echo "0"
echo "11"
if ($s2t == 1) then
	ev "$s1c*($rg_r_arr[$s1t]-$rg_r_arr[1])" \
		"$s1c*($rg_g_arr[$s1t]-$rg_g_arr[1])" \
		"$s1c*($rg_b_arr[$s1t]-$rg_b_arr[1])"
	ev "$s1c*$rc_r_arr[$s1t]" "$s1c*$rc_g_arr[$s1t]" "$s1c*$rc_b_arr[$s1t]"
	ev "$s1c*$tn_r_arr[$s1t]" "$s1c*$tn_g_arr[$s1t]" "$s1c*$tn_b_arr[$s1t]"
	echo "	1	$s1c"
else
	ev "$s2c*$rc_r_arr[$s2t]" "$s2c*$rc_g_arr[$s2t]" "$s2c*$rc_b_arr[$s2t]"
	ev "$s2c*($rg_r_arr[$s2t]-$rg_r_arr[1])" \
		"$s2c*($rg_g_arr[$s2t]-$rg_g_arr[1])" \
		"$s2c*($rg_b_arr[$s2t]-$rg_b_arr[1])"
	ev "$s2c*$tn_r_arr[$s2t]" "$s2c*$tn_g_arr[$s2t]" "$s2c*$tn_b_arr[$s2t]"
	echo "	-1	$s2c"
endif
else
### Low-E glazing
echo "10"
echo "	sr_lowE_r sr_lowE_g sr_lowE_b"
echo "	st_lowE_r st_lowE_g st_lowE_b"
echo "	0 0 0"
echo "	glaze1.cal"
echo "0"
echo "19"
echo "	0 0 0"
echo "	0 0 0"
echo "	0 0 0"
if ($s2t == 1) then
	echo "	1"
	set st=$s1t
else
	echo "	-1"
	set st=$s2t
endif
echo "	$rg_r_arr[$st] $rg_g_arr[$st] $rg_b_arr[$st]"
echo "	$rc_r_arr[$st] $rc_g_arr[$st] $rc_b_arr[$st]"
echo "	$tn_r_arr[$st] $tn_g_arr[$st] $tn_b_arr[$st]"
endif
echo ""
exit 0

glaze2:
#################################################################
#
# Compute double glazing
#
if ($s2t != 1) then
	set s2r_rgb=($rc_r_arr[$s2t] $rc_g_arr[$s2t] $rc_b_arr[$s2t])
	set s1r_rgb=($rg_r_arr[$s2t] $rg_g_arr[$s2t] $rg_b_arr[$s2t])
	set s12t_rgb=($tn_r_arr[$s2t] $tn_g_arr[$s2t] $tn_b_arr[$s2t])
else
	set s2r_rgb=($rg_r_arr[$s1t] $rg_g_arr[$s1t] $rg_b_arr[$s1t])
	set s1r_rgb=($rc_r_arr[$s1t] $rc_g_arr[$s1t] $rc_b_arr[$s1t])
	set s12t_rgb=($tn_r_arr[$s1t] $tn_g_arr[$s1t] $tn_b_arr[$s1t])
endif
if ($s4t != 1) then
	set s4r_rgb=($rc_r_arr[$s4t] $rc_g_arr[$s4t] $rc_b_arr[$s4t])
	set s3r_rgb=($rg_r_arr[$s4t] $rg_g_arr[$s4t] $rg_b_arr[$s4t])
	set s34t_rgb=($tn_r_arr[$s4t] $tn_g_arr[$s4t] $tn_b_arr[$s4t])
else
	set s4r_rgb=($rg_r_arr[$s3t] $rg_g_arr[$s3t] $rg_b_arr[$s3t])
	set s3r_rgb=($rc_r_arr[$s3t] $rc_g_arr[$s3t] $rc_b_arr[$s3t])
	set s34t_rgb=($tn_r_arr[$s3t] $tn_g_arr[$s3t] $tn_b_arr[$s3t])
endif
set s12c=1
echo "# Exterior surface s1 type: $sn_arr[$s1t]"
if ($?s1c) then
	echo "# s1 coating coverage: $s1c"
	set s12c=$s1c
endif
echo "# Inner surface s2 type: $sn_arr[$s2t]"
if ($?s2c) then
	echo "# s2 coating coverage: $s2c"
	set s12c=$s2c
endif
set s34c=1
echo "# Inner surface s3 type: $sn_arr[$s3t]"
if ($?s3c) then
	echo "# s3 coating coverage: $s3c"
	set s34c=$s3c
endif
echo "# Interior surface s4 type: $sn_arr[$s4t]"
if ($?s4c) then
	echo "# s4 coating coverage: $s4c"
	set s34c=$s4c
endif
# Approximate reflectance and transmittance for comment using gray values
set rglass=`ev ".265*$rg_r_arr[1]+.670*$rg_g_arr[1]+.065*$rg_b_arr[1]"`
set tglass=`ev ".265*$tn_r_arr[1]+.670*$tn_g_arr[1]+.065*$tn_b_arr[1]"`
set s1r_gry=`ev "$s12c*(.265*$s1r_rgb[1]+.670*$s1r_rgb[2]+.065*$s1r_rgb[3])+(1-$s12c)*$rglass"`
set s2r_gry=`ev "$s12c*(.265*$s2r_rgb[1]+.670*$s2r_rgb[2]+.065*$s2r_rgb[3])+(1-$s12c)*$rglass"`
set s12t_gry=`ev "$s12c*(.265*$s12t_rgb[1]+.670*$s12t_rgb[2]+.065*$s12t_rgb[3])+(1-$s12c)*$tglass"`
set s3r_gry=`ev "$s34c*(.265*$s3r_rgb[1]+.670*$s3r_rgb[2]+.065*$s3r_rgb[3])+(1-$s34c)*$rglass"`
set s4r_gry=`ev "$s34c*(.265*$s4r_rgb[1]+.670*$s4r_rgb[2]+.065*$s4r_rgb[3])+(1-$s34c)*$rglass"`
set s34t_gry=`ev "$s34c*(.265*$s34t_rgb[1]+.670*$s34t_rgb[2]+.065*$s34t_rgb[3])+(1-$s34c)*$tglass"`
echo -n "# Exterior normal hemispherical reflectance: "
ev "$s1r_gry + $s12t_gry^2*$s3r_gry"
echo -n "# Interior normal hemispherical reflectance: "
ev "$s4r_gry + $s34t_gry^2*$s2r_gry"
echo -n "# Normal hemispherical transmittance: "
ev "$s12t_gry*$s34t_gry"
echo "#"
echo "void BRTDfunc glaze2_unnamed"

if ($part_arr[$s3t] || $part_arr[$s4t]) then
### Front pane has frit
if ($part_arr[$s1t] || $part_arr[$s2t]) then
	echo "Only one pane can have frit"
	exit 1
endif
if ($?s3c) then
	set sc=$s3c
	set s3g=`ev "1-$s3c"`
else
	set s3c=0
	set s3g=1
endif
if ($?s4c) then
	set sc=$s4c
	set s4g=`ev "1-$s4c"`
else
	set s4c=0
	set s4g=1
endif
echo "10"
echo "if(Rdot,cr($s4g*rclr,$s3g*$s4g*tclr,fr($s2r_rgb[1])),cr(fr($s1r_rgb[1]),ft($s12t_rgb[1]),$s3g*rclr))"
echo "if(Rdot,cr($s4g*rclr,$s3g*$s4g*tclr,fr($s2r_rgb[2])),cr(fr($s1r_rgb[2]),ft($s12t_rgb[2]),$s3g*rclr))"
echo "if(Rdot,cr($s4g*rclr,$s3g*$s4g*tclr,fr($s2r_rgb[3])),cr(fr($s1r_rgb[3]),ft($s12t_rgb[3]),$s3g*rclr))"
echo "$s3g*$s4g*ft($s12t_rgb[1])*tclr"
echo "$s3g*$s4g*ft($s12t_rgb[2])*tclr"
echo "$s3g*$s4g*ft($s12t_rgb[3])*tclr"
echo "	0 0 0"
echo "	glaze2.cal"
echo "0"
echo "9"
ev "$sc*$s4r_rgb[1]-$s3c*$rg_r_arr[1]" \
	"$sc*$s4r_rgb[2]-$s3c*$rg_g_arr[1]" \
	"$sc*$s4r_rgb[3]-$s3c*$rg_b_arr[1]"
ev "$s12t_rgb[1]^2*($sc*$s3r_rgb[1]-$s4c*$rg_r_arr[1])" \
	"$s12t_rgb[2]^2*($sc*$s3r_rgb[2]-$s4c*$rg_g_arr[1])" \
	"$s12t_rgb[3]^2*($sc*$s3r_rgb[3]-$s4c*$rg_b_arr[1])"
ev "$sc*$s12t_rgb[1]*$s34t_rgb[1]" \
	"$sc*$s12t_rgb[2]*$s34t_rgb[2]" \
	"$sc*$s12t_rgb[3]*$s34t_rgb[3]"
else if ($part_arr[$s1t] || $part_arr[$s2t]) then
### Back pane has frit
if ($?s1c) then
	set sc=$s1c
	set s1g=`ev "1-$s1c"`
else
	set s1c=0
	set s1g=1
endif
if ($?s2c) then
	set sc=$s2c
	set s2g=`ev "1-$s2c"`
else
	set s2c=0
	set s2g=1
endif
echo "10"
echo "if(Rdot,cr(fr($s4r_rgb[1]),ft($s34t_rgb[1]),$s2g*rclr),cr($s1g*rclr,$s1g*$s2g*tclr,fr($s3r_rgb[1])))"
echo "if(Rdot,cr(fr($s4r_rgb[2]),ft($s34t_rgb[2]),$s2g*rclr),cr($s1g*rclr,$s1g*$s2g*tclr,fr($s3r_rgb[2])))"
echo "if(Rdot,cr(fr($s4r_rgb[3]),ft($s34t_rgb[3]),$s2g*rclr),cr($s1g*rclr,$s1g*$s2g*tclr,fr($s3r_rgb[3])))"
echo "$s1g*$s2g*ft($s34t_rgb[1])*tclr"
echo "$s1g*$s2g*ft($s34t_rgb[2])*tclr"
echo "$s1g*$s2g*ft($s34t_rgb[3])*tclr"
echo "	0 0 0"
echo "	glaze2.cal"
echo "0"
echo "9"
ev "$s34t_rgb[1]^2*($sc*$s2r_rgb[1]-$s1c*$rg_r_arr[1])" \
	"$s34t_rgb[2]^2*($sc*$s2r_rgb[2]-$s1c*$rg_g_arr[1])" \
	"$s34t_rgb[3]^2*($sc*$s2r_rgb[3]-$s1c*$rg_b_arr[1])"
ev "$sc*$s1r_rgb[1]-$s2c*$rg_r_arr[1]" \
	"$sc*$s1r_rgb[2]-$s2c*$rg_g_arr[1]" \
	"$sc*$s1r_rgb[3]-$s2c*$rg_b_arr[1]"
ev "$sc*$s34t_rgb[1]*$s12t_rgb[1]" \
	"$sc*$s34t_rgb[2]*$s12t_rgb[2]" \
	"$sc*$s34t_rgb[3]*$s12t_rgb[3]"
else
### Low-E and regular glazing only
echo "10"
echo "if(Rdot,cr(fr($s4r_rgb[1]),ft($s34t_rgb[1]),fr($s2r_rgb[1])),cr(fr($s1r_rgb[1]),ft($s12t_rgb[1]),fr($s3r_rgb[1])))"
echo "if(Rdot,cr(fr($s4r_rgb[2]),ft($s34t_rgb[2]),fr($s2r_rgb[2])),cr(fr($s1r_rgb[2]),ft($s12t_rgb[2]),fr($s3r_rgb[2])))"
echo "if(Rdot,cr(fr($s4r_rgb[3]),ft($s34t_rgb[3]),fr($s2r_rgb[3])),cr(fr($s1r_rgb[3]),ft($s12t_rgb[3]),fr($s3r_rgb[3])))"
echo "ft($s34t_rgb[1])*ft($s12t_rgb[1])"
echo "ft($s34t_rgb[2])*ft($s12t_rgb[2])"
echo "ft($s34t_rgb[3])*ft($s12t_rgb[3])"
echo "	0 0 0"
echo "	glaze2.cal"
echo "0"
echo "9"
echo "	0 0 0"
echo "	0 0 0"
echo "	0 0 0"
endif
echo ""
exit 0
