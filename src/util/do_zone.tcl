# RCSid: $Id: do_zone.tcl,v 2.7 2003/02/22 02:07:30 greg Exp $
#
# Create ZONE screen
#

proc zonechange {nm el op} {		# record change in ZONE
	global radvar zonevar
	if {"$nm" == "radvar"} {
		if {"$op" == "u"} {
			set radvar(ZONE) {}
			trace variable radvar(ZONE) wu zonechange
		}
		trace vdelete zonevar w zonechange
		if {"$radvar(ZONE)" == {}} {
			set zonevar(IE) Interior
			set zonevar(xmin) {}
			set zonevar(xmax) {}
			set zonevar(ymin) {}
			set zonevar(ymax) {}
			set zonevar(zmin) {}
			set zonevar(zmax) {}
		} elseif {[llength $radvar(ZONE)] == 7} {
			if [string match {[Ii]*} $radvar(ZONE)] {
				set zonevar(IE) Interior
			} else {
				set zonevar(IE) Exterior
			}
			set zonevar(xmin) [lindex $radvar(ZONE) 1]
			set zonevar(xmax) [lindex $radvar(ZONE) 2]
			set zonevar(ymin) [lindex $radvar(ZONE) 3]
			set zonevar(ymax) [lindex $radvar(ZONE) 4]
			set zonevar(zmin) [lindex $radvar(ZONE) 5]
			set zonevar(zmax) [lindex $radvar(ZONE) 6]
		}
		trace variable zonevar w zonechange
		return
	}
	trace vdelete radvar(ZONE) wu zonechange
	set radvar(ZONE) "$zonevar(IE) $zonevar(xmin) $zonevar(xmax) $zonevar(ymin) $zonevar(ymax) $zonevar(zmin) $zonevar(zmax)"
	trace variable radvar(ZONE) wu zonechange
}

proc indscale v {		# set indirect value
	global radvar
	set radvar(INDIRECT) $v
}

proc copyzone rf {		# copy ZONE settings from another RIF
	load_vars [file tail $rf] {ZONE DETAIL INDIRECT VARIABILITY EXPOSURE}
}

proc add2bb f {
	global bbargs
	lappend bbargs $f
}

proc get_bbox {} {
	global radvar zonevar myglob bbargs curmess curpat
	set bbargs {}
	if [getfile -grab -perm -glob $myglob(scene) -view view_txt -send add2bb] {
		set curmess "Running getbbox..."
		update idletasks
		if [catch {eval exec getbbox -h $bbargs} bbout] {
			set curmess $bbout
		} else {
			set curmess "Done."
			set radvar(ZONE) "$zonevar(IE) $bbout"
		}
	}
	set myglob(scene) $curpat
}

proc do_zone w {		# set up ZONE screen
	global radvar zonevar rifname
	if {"$w" == "done"} {
		trace vdelete radvar(ZONE) wu zonechange
		trace vdelete zonevar w zonechange
		unset zonevar
		return
	}
	frame $w
	# Zone fields
	label $w.zol -text ZONE
	place $w.zol -relx .0714 -rely .0610
	radiobutton $w.zi -text Interior -anchor w -relief flat \
			-variable zonevar(IE) -value Interior
	place $w.zi -relx .2857 -rely .0610
	radiobutton $w.ze -text Exterior -anchor w -relief flat \
			-variable zonevar(IE) -value Exterior
	place $w.ze -relx .5714 -rely .0610
	helplink "$w.zi $w.ze" trad zone type
	label $w.maxl -text Maximum
	place $w.maxl -relx .0714 -rely .1829
	label $w.minl -text Minimum
	place $w.minl -relx .0714 -rely .3049
	button $w.bbox -text Auto -relief raised -command get_bbox
	place $w.bbox -relwidth .1071 -relheight .0610 -relx .0714 -rely .237
	helplink $w.bbox trad zone auto
	label $w.xl -text X:
	place $w.xl -relx .2500 -rely .2439
	entry $w.xmax -textvariable zonevar(xmax) -relief sunken
	place $w.xmax -relwidth .1786 -relheight .0610 -relx .2857 -rely .1829
	entry $w.xmin -textvariable zonevar(xmin) -relief sunken
	place $w.xmin -relwidth .1786 -relheight .0610 -relx .2857 -rely .3049
	label $w.yl -text Y:
	place $w.yl -relx .5000 -rely .2439
	entry $w.ymax -textvariable zonevar(ymax) -relief sunken
	place $w.ymax -relwidth .1786 -relheight .0610 -relx .5357 -rely .1829
	entry $w.ymin -textvariable zonevar(ymin) -relief sunken
	place $w.ymin -relwidth .1786 -relheight .0610 -relx .5357 -rely .3049
	label $w.zl -text Z:
	place $w.zl -relx .7500 -rely .2439
	entry $w.zmax -textvariable zonevar(zmax) -relief sunken
	place $w.zmax -relwidth .1786 -relheight .0610 -relx .7857 -rely .1829
	entry $w.zmin -textvariable zonevar(zmin) -relief sunken
	place $w.zmin -relwidth .1786 -relheight .0610 -relx .7857 -rely .3049
	helplink "$w.xmin $w.xmax $w.ymin $w.ymax $w.zmin $w.zmax" \
			trad zone zone
	# Detail
	label $w.detl -text Detail:
	place $w.detl -relx .0714 -rely .4268
	radiobutton $w.detH -text High -anchor w -relief flat \
			-variable radvar(DETAIL) -value High
	place $w.detH -relx .2857 -rely .4268
	radiobutton $w.detM -text Medium -anchor w -relief flat \
			-variable radvar(DETAIL) -value Medium
	place $w.detM -relx .5357 -rely .4268
	radiobutton $w.detL -text Low -anchor w -relief flat \
			-variable radvar(DETAIL) -value Low
	place $w.detL -relx .7857 -rely .4268
	helplink "$w.detH $w.detM $w.detL" trad zone detail
	# Indirect
	label $w.indl -text Indirect:
	place $w.indl -relx .0714 -rely .5488
	scale $w.inds -showvalue yes -from 0 -to 5 \
			-orient horizontal -command indscale
	$w.inds set $radvar(INDIRECT)
	place $w.inds -relx .2857 -rely .5488 -relwidth .6500 -relheight .1200
	helplink $w.inds trad zone indirect
	# Variability
	label $w.varl -text Variability:
	place $w.varl -relx .0714 -rely .7317
	radiobutton $w.varH -text High -anchor w -relief flat \
			-variable radvar(VARIABILITY) -value High
	place $w.varH -relx .2857 -rely .7317
	radiobutton $w.varM -text Medium -anchor w -relief flat \
			-variable radvar(VARIABILITY) -value Medium
	place $w.varM -relx .5357 -rely .7317
	radiobutton $w.varL -text Low -anchor w -relief flat \
			-variable radvar(VARIABILITY) -value Low
	place $w.varL -relx .7857 -rely .7317
	helplink "$w.varH $w.varM $w.varL" trad zone variability
	# Exposure
	label $w.expl -text Exposure:
	place $w.expl -relx .0714 -rely .8537
	entry $w.expe -textvariable radvar(EXPOSURE) -relief sunken
	place $w.expe -relwidth .1786 -relheight .0610 -relx .2857 -rely .8537
	helplink $w.expe trad zone exposure
	# Revert and Copy buttons
	button $w.revert -text Revert -relief raised \
			-command "copyzone $rifname"
	place $w.revert -relwidth .1071 -relheight .0610 -relx .98 -rely .98 \
			-anchor se
	helplink $w.revert trad zone revert
	button $w.copy -text Copy -relief raised -command {getfile -grab \
			-send copyzone -view view_txt -glob $myglob(rif)}
	place $w.copy -relwidth .1071 -relheight .0610 -relx .98 -rely .90 \
			-anchor se
	helplink $w.copy trad zone copy
	# Focusing
	focus $w.xmin
	bind $w.xmin <Return> "focus $w.xmax"
	bind $w.xmax <Return> "focus $w.ymin"
	bind $w.ymin <Return> "focus $w.ymax"
	bind $w.ymax <Return> "focus $w.zmin"
	bind $w.zmin <Return> "focus $w.zmax"
	bind $w.zmax <Return> "focus $w.expe"
	bind $w.expe <Return> {}
	# Track ZONE variable
	if [info exists radvar(ZONE)] {
		zonechange radvar ZONE w
		trace variable radvar(ZONE) wu zonechange
	} else {
		zonechange radvar ZONE u
	}
}
