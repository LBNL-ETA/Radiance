# RCSid: $Id: do_views3.6.tcl,v 2.8 2003/02/22 02:07:30 greg Exp $
#
# Handle view parameters
#

proc viewchange {nm el op} {		# handle changes to radvar(view)
	global radvar viewbox viewname viewopts
	if {"$op" == "u"} {
		set radvar(view) {}
		trace variable radvar(view) wu viewchange
	}
	$viewbox delete 0 end
	foreach v $radvar(view) {
		$viewbox insert end [lindex $v 0]
	}
}

proc vnchange {nm el op} {		# handle changes to viewname
	global viewname vdefbutt radvar
	if {"$viewname" == {}} {
		$vdefbutt configure -text "Set Default" -state disabled
	} elseif {$radvar(view) != {} &&
			"[lindex [lindex $radvar(view) 0] 0]" == "$viewname"} {
		$vdefbutt configure -text "Is Default" -state disabled
	} else {
		$vdefbutt configure -text "Set Default" -state normal
	}
}

proc selview {} {			# change selected view
	global radvar viewbox viewname viewopts viewseln
	set viewseln [$viewbox curselection]
	set viewname [$viewbox get $viewseln]
	set viewopts [lrange [lindex $radvar(view) $viewseln] 1 end]
}

proc addview {{pos end}} {		# add current view
	global radvar viewname viewopts curmess viewseln viewbox
	if {$pos == -1} {return -1}
	if {[llength $viewname] != 1} {
		beep
		set curmess "Illegal view name."
		return -1
	}
	set curmess {}
	set n [lsearch -exact $radvar(view) $viewname]
	if {$n < 0} {
		set n [lsearch -glob $radvar(view) "$viewname *"]
	}
	if {$n >= 0} {
		set radvar(view) [lreplace $radvar(view) $n $n]
		if {"$pos" != "end" && $pos > $n} {incr pos -1}
	}
	if {"$pos" == "end"} {
		set viewseln [llength $radvar(view)]
		lappend radvar(view) "$viewname $viewopts"
	} else {
		set viewseln $pos
		set radvar(view) [linsert $radvar(view) $pos \
				"$viewname $viewopts"]
	}
	vnchange viewname {} w
	$viewbox select from $pos
	return $pos
}

proc delview {} {			# delete current view
	global radvar viewbox curmess viewseln
	if {! [info exists viewseln]} {
		beep
		set curmess "No view selected for deletion."
		return -1
	}
	set n $viewseln
	set curmess {}
	set radvar(view) [lreplace $radvar(view) $n $n]
	vnchange viewname {} w
	unset viewseln
	return $n
}

proc copyviews rf {		# copy view settings from another RIF
	load_vars [file tail $rf] {view UP PICTURE RESOLUTION RAWFILE ZFILE}
	vnchange viewname {} w
}

proc do_views w {			# create views screen
	global radvar viewbox viewseln viewname viewopts vdefbutt rifname
	if {"$w" == "done"} {
		trace vdelete radvar(view) wu viewchange
		trace vdelete viewname w vnchange
		unset viewbox viewname viewopts vdefbutt
		return
	}
	frame $w
	# View box
	label $w.vvl -text Views
	place $w.vvl -relx .0714 -rely .0610
	frame $w.vnl
	scrollbar $w.vnl.sb -relief sunken -command "$w.vnl.lb yview"
	listbox $w.vnl.lb -relief sunken -yscroll "$w.vnl.sb set"
	set viewbox $w.vnl.lb
	tk_listboxSingleSelect $w.vnl.lb
	bind $w.vnl.lb <ButtonRelease-1> +selview
	pack $w.vnl.sb -side right -fill y
	pack $w.vnl.lb -side left -expand yes -fill both
	place $w.vnl -relwidth .1786 -relheight .3049 -relx .07143 -rely .1220
	helplink $w.vnl.lb trad views list
	# View name
	set viewname {}
	label $w.vnlb -text Name
	place $w.vnlb -relx .2857 -rely .0610
	entry $w.vne -textvariable viewname -relief sunken
	place $w.vne -relwidth .1429 -relheight .0610 -relx .2857 -rely .1220
	helplink $w.vne trad views name
	# View options
	set viewopts {}
	label $w.vol -text Options
	place $w.vol -relx .2857 -rely .1829
	frame $w.vo
	entry $w.vo.e -textvariable viewopts -relief sunken \
			-scroll "$w.vo.s set"
	scrollbar $w.vo.s -relief sunken -orient horiz -command "$w.vo.e view"
	pack $w.vo.e $w.vo.s -side top -fill x
	place $w.vo -relwidth .6429 -relheight .1000 -relx .2857 -rely .2439
	helplink $w.vo.e trad views options
	# View change buttons
	button $w.vclr -text Clear -width 7 -relief raised \
			-command "set viewname {} ; set viewopts {}
					focus $w.vne
					catch {unset viewseln}
					$w.vnl.lb select clear"
	place $w.vclr -relx .4643 -rely .1220
	helplink $w.vclr trad views clear
	button $w.vdef -text "Set Default" -width 12 \
			-relief raised -state disabled -command {addview 0}
	place $w.vdef -relx .6429 -rely .1220
	helplink $w.vdef trad views default
	set vdefbutt $w.vdef
	button $w.vadd -text Add -width 7 -relief raised -command addview
	place $w.vadd -relx .2857 -rely .3659
	helplink $w.vadd trad views add
	button $w.vchg -text Change -width 7 -relief raised \
			-command {addview [delview]}
	place $w.vchg -relx .4643 -rely .3659
	helplink $w.vchg trad views change
	button $w.vdel -text Delete -width 7 -relief raised -command delview
	place $w.vdel -relx .6429 -rely .3659
	helplink $w.vdel trad views delete
	# Up vector
	label $w.upl -text Up:
	place $w.upl -relx .0714 -rely .4878
	radiobutton $w.upxp -text +X -variable radvar(UP) -value X -relief flat
	place $w.upxp -relx .2857 -rely .4634
	radiobutton $w.upxn -text -X -variable radvar(UP) -value -X -relief flat
	place $w.upxn -relx .2857 -rely .5244
	radiobutton $w.upyp -text +Y -variable radvar(UP) -value Y -relief flat
	place $w.upyp -relx .4286 -rely .4634
	radiobutton $w.upyn -text -Y -variable radvar(UP) -value -Y -relief flat
	place $w.upyn -relx .4286 -rely .5244
	radiobutton $w.upzp -text +Z -variable radvar(UP) -value Z -relief flat
	place $w.upzp -relx .5714 -rely .4634
	radiobutton $w.upzn -text -Z -variable radvar(UP) -value -Z -relief flat
	place $w.upzn -relx .5714 -rely .5244
	helplink "$w.upxp $w.upxn $w.upyp $w.upyn $w.upzp $w.upzn" \
			trad views up
	# Picture file name
	label $w.pfl -text Picture:
	place $w.pfl -relx .0714 -rely .6345
	entry $w.pfe -textvariable radvar(PICTURE) -relief sunken
	place $w.pfe -relwidth .5714 -relheight .0610 -relx .2857 -rely .6345
	helplink $w.pfe trad views picture
	# Resolution
	label $w.rzl -text Resolution:
	place $w.rzl -relx .0714 -rely .7107
	entry $w.rze -textvariable radvar(RESOLUTION) -relief sunken
	place $w.rze -relwidth .2857 -relheight .0610 -relx .2857 -rely .7107
	helplink $w.rze trad views resolution
	# Rawfile
	label $w.rfl -text Rawfile:
	place $w.rfl -relx .0714 -rely .7870
	entry $w.rfe -textvariable radvar(RAWFILE) -relief sunken
	place $w.rfe -relwidth .5714 -relheight .0610 -relx .2857 -rely .7870
	helplink $w.rfe trad views rawfile
	# Zfile
	label $w.zfl -text Zfile:
	place $w.zfl -relx .0714 -rely .8632
	entry $w.zfe -textvariable radvar(ZFILE) -relief sunken
	place $w.zfe -relwidth .5714 -relheight .0610 -relx .2857 -rely .8632
	helplink $w.zfe trad views zfile
	# Revert and Copy buttons
	button $w.revert -text Revert -relief raised \
			-command "copyviews $rifname"
	place $w.revert -relwidth .1071 -relheight .0610 -relx .98 -rely .98 \
			-anchor se
	helplink $w.revert trad views revert
	button $w.copy -text Copy -relief raised -command {getfile -grab \
			-send copyviews -view view_txt -glob $myglob(rif)}
	place $w.copy -relwidth .1071 -relheight .0610 -relx .98 -rely .90 \
			-anchor se
	helplink $w.copy trad views copy
	# Assign focus
	bind $w.vne <Return> "focus $w.vo.e"
	bind $w.vo.e <Return> "addview ; $w.vclr invoke"
	bind $w.vo.e <B2-Motion> {}
	bind $w.vo.e <Button-2> "$w.vo.e insert insert \[selection get\]"
	$w.vclr invoke
	# Track variables
	vnchange viewname {} w
	trace variable viewname w vnchange
	if [info exists radvar(view)] {
		viewchange radvar view w
		trace variable radvar(view) wu viewchange
	} else {
		viewchange radvar view u
	}
}
