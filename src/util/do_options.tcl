# RCSid: $Id: do_options.tcl,v 2.8 2003/02/22 02:07:30 greg Exp $
#
# Options screen for trad
#

proc amb_delete {} {		# delete ambient file
	global radvar curmess
	if {"$radvar(AMBFILE)" == {} || ! [file isfile $radvar(AMBFILE)]} {
		set curmess {No ambient file.}
		return
	}
	set hl [string length [exec getinfo < $radvar(AMBFILE)]]
	set fl [file size $radvar(AMBFILE)]
	set nv [expr ($fl - $hl - 2) / 75]
	if {$nv > 50 && [tk_dialog .dlg {Verification} \
				"Really delete ambient file $radvar(AMBFILE)\
				with $nv indirect irradiance values?" \
				questhead 0 {Delete} {Cancel}]} {
		return 0
	}
	if [catch {exec rm $radvar(AMBFILE) < /dev/null} curmess] {return 0}
	set curmess {Ambient file deleted.}
	return 1
}

proc copyopts rf {		# copy option settings from another RIF
	load_vars [file tail $rf] {QUALITY PENUMBRAS AMBFILE OPTFILE REPORT \
					oconv mkillum render pfilt}
}

proc do_options w {		# Set up options screen
	global radvar rifname
	if {"$w" == "done"} {
		return
	}
	frame $w
	# Quality
	label $w.qual -text Quality:
	place $w.qual -relx .0714 -rely .0610
	radiobutton $w.quaH -text High -anchor w -relief flat \
			-variable radvar(QUALITY) -value High
	place $w.quaH -relx .2857 -rely .0610
	radiobutton $w.quaM -text Medium -anchor w -relief flat \
			-variable radvar(QUALITY) -value Medium
	place $w.quaM -relx .5357 -rely .0610
	radiobutton $w.quaL -text Low -anchor w -relief flat \
			-variable radvar(QUALITY) -value Low
	place $w.quaL -relx .7857 -rely .0610
	helplink "$w.quaH $w.quaM $w.quaL" trad options quality
	# Penumbras
	label $w.penl -text Penumbras:
	place $w.penl -relx .0714 -rely .1585
	radiobutton $w.penT -text On -anchor w -relief flat \
			-variable radvar(PENUMBRAS) -value True
	place $w.penT -relx .2857 -rely .1585
	radiobutton $w.penF -text Off -anchor w -relief flat \
			-variable radvar(PENUMBRAS) -value False
	place $w.penF -relx .5357 -rely .1585
	helplink "$w.penT $w.penF" trad options penumbras
	# Ambfile
	label $w.ambl -text Ambfile:
	place $w.ambl -relx .0714 -rely .2805
	entry $w.ambe -textvariable radvar(AMBFILE) -relief sunken
	place $w.ambe -relwidth .5714 -relheight .0610 -relx .2857 -rely .2805
	helplink $w.ambe trad options ambfile
	button $w.ambd -text Delete -relief raised -command amb_delete
	place $w.ambd -relwidth .1071 -relheight .0610 -relx .98 -rely .2805 \
			-anchor ne
	helplink $w.ambd trad options ambdelete
	# Optfile
	label $w.optl -text Optfile:
	place $w.optl -relx .0714 -rely .4024
	entry $w.opte -textvariable radvar(OPTFILE) -relief sunken
	place $w.opte -relwidth .5714 -relheight .0610 -relx .2857 -rely .4024
	helplink $w.opte trad options optfile
	# Report
	label $w.repl -text Report:
	place $w.repl -relx .0714 -rely .5244
	entry $w.repe -textvariable radvar(REPORT) -relief sunken
	place $w.repe -relwidth .5714 -relheight .0610 -relx .2857 -rely .5244
	helplink $w.repe trad options report
	# Oconv options
	label $w.ocol -text "oconv opts:"
	place $w.ocol -relx .0714 -rely .6463
	entry $w.ocoe -textvariable radvar(oconv) -relief sunken
	place $w.ocoe -relwidth .5714 -relheight .0600 -relx .2857 -rely .6463
	helplink $w.ocoe trad options oconv
	# Mkillum options
	label $w.mkil -text "mkillum opts:"
	place $w.mkil -relx .0714 -rely .7073
	entry $w.mkie -textvariable radvar(mkillum) -relief sunken
	place $w.mkie -relwidth .5714 -relheight .0600 -relx .2857 -rely .7073
	helplink $w.mkie trad options mkillum
	# Render options
	label $w.renl -text "render opts:"
	place $w.renl -relx .0714 -rely .7683
	entry $w.rene -textvariable radvar(render) -relief sunken
	place $w.rene -relwidth .5714 -relheight .0600 -relx .2857 -rely .7683
	helplink $w.rene trad options render
	# Pfilt options
	label $w.pfil -text "pfilt opts:"
	place $w.pfil -relx .0714 -rely .8293
	entry $w.pfie -textvariable radvar(pfilt) -relief sunken
	place $w.pfie -relwidth .5714 -relheight .0600 -relx .2857 -rely .8293
	helplink $w.pfie trad options pfilt
	# Revert and Copy buttons
	button $w.revert -text Revert -relief raised \
			-command "copyopts $rifname"
	place $w.revert -relwidth .1071 -relheight .0610 -relx .98 -rely .98 \
			-anchor se
	helplink $w.revert trad options revert
	button $w.copy -text Copy -relief raised -command {getfile -grab \
			-send copyopts -view view_txt -glob $myglob(rif)}
	place $w.copy -relwidth .1071 -relheight .0610 -relx .98 -rely .90 \
			-anchor se
	helplink $w.copy trad options copy
	# Set focus
	focus $w.ambe
	bind $w.ambe <Return> "focus $w.opte"
	bind $w.opte <Return> "focus $w.repe"
	bind $w.repe <Return> "focus $w.ocoe"
	bind $w.ocoe <Return> "focus $w.mkie"
	bind $w.mkie <Return> "focus $w.rene"
	bind $w.rene <Return> "focus $w.pfie"
	bind $w.pfie <Return> {}
}
