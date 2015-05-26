# RCSid: $Id: do_options.tcl,v 2.10 2015/05/26 12:39:25 greg Exp $
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
	set nv [expr ($fl - $hl - 2) / 67]
	if {$nv > 1000 && [tk_dialog .dlg {Verification} \
				"Really delete ambient file $radvar(AMBFILE)\
				with $nv indirect irradiance values?" \
				questhead 0 {Delete} {Cancel}]} {
		return 0
	}
	if [catch {exec rm $radvar(AMBFILE) < /dev/null} curmess] {return 0}
	set curmess {Ambient file deleted.}
	return 1
}

proc pgm_delete {} {		# delete global photon map file
	global radvar curmess
	set pm [lindex $radvar(PGMAP) 0]
	if {"$pm" == {} || ! [file isfile $pm]} {
		set curmess {No global photon file.}
		return
	}
	set hd [exec getinfo < $pm]
	set nv [lindex $hd [expr [lsearch $hd photons] - 1]]
	if {$nv > 10000 && [tk_dialog .dlg {Verification} \
				"Really delete global photon map $pm\
				with $nv photons?" \
				questhead 0 {Delete} {Cancel}]} {
		return 0
	}
	if [catch {exec rm $pm < /dev/null} curmess] {return 0}
	set curmess {Global photon file deleted.}
	return 1
}

proc pcm_delete {} {		# delete caustic photon map file
	global radvar curmess
	set pm [lindex $radvar(PCMAP) 0]
	if {"$pm" == {} || ! [file isfile $pm]} {
		set curmess {No caustic photon file.}
		return
	}
	set hd [exec getinfo < $pm]
	set nv [lindex $hd [expr [lsearch $hd photons] - 1]]
	if {$nv > 10000 && [tk_dialog .dlg {Verification} \
				"Really delete caustic photon map $pm\
				with $nv photons?" \
				questhead 0 {Delete} {Cancel}]} {
		return 0
	}
	if [catch {exec rm $pm < /dev/null} curmess] {return 0}
	set curmess {Caustic photon file deleted.}
	return 1
}

proc copyopts rf {		# copy option settings from another RIF
	load_vars [file tail $rf] {QUALITY PENUMBRAS AMBFILE PGMAP PCMAP OPTFILE REPORT \
					oconv mkillum mkpmap render pfilt}
}

proc do_options w {		# Set up options screen
	global radvar rifname
	if {"$w" == "done"} {
		return
	}
	frame $w
	# Quality
	label $w.qual -text Quality:
	place $w.qual -relx .0714 -rely .0510
	radiobutton $w.quaH -text High -anchor w -relief flat \
			-variable radvar(QUALITY) -value High
	place $w.quaH -relx .2857 -rely .0510
	radiobutton $w.quaM -text Medium -anchor w -relief flat \
			-variable radvar(QUALITY) -value Medium
	place $w.quaM -relx .5357 -rely .0510
	radiobutton $w.quaL -text Low -anchor w -relief flat \
			-variable radvar(QUALITY) -value Low
	place $w.quaL -relx .7857 -rely .0510
	helplink "$w.quaH $w.quaM $w.quaL" trad options quality
	# Penumbras
	label $w.penl -text Penumbras:
	place $w.penl -relx .0714 -rely .1485
	radiobutton $w.penT -text On -anchor w -relief flat \
			-variable radvar(PENUMBRAS) -value True
	place $w.penT -relx .2857 -rely .1485
	radiobutton $w.penF -text Off -anchor w -relief flat \
			-variable radvar(PENUMBRAS) -value False
	place $w.penF -relx .5357 -rely .1485
	helplink "$w.penT $w.penF" trad options penumbras
	# Global photon map
	label $w.pgml -text Pgmap:
	place $w.pgml -relx .0714 -rely .2460
	entry $w.pgme -textvariable radvar(PGMAP) -relief sunken
	place $w.pgme -relwidth .5714 -relheight .0610 -relx .2857 -rely .2460
	helplink $w.pgme trad options pgmap
	button $w.pgmd -text Delete -relief raised -command pgm_delete
	place $w.pgmd -relwidth .1071 -relheight .0610 -relx .98 -rely .2460 \
			-anchor ne
	helplink $w.pgmd trad options pgmdelete
	# Caustic photon map
	label $w.pcml -text Pcmap:
	place $w.pcml -relx .0714 -rely .3341
	entry $w.pcme -textvariable radvar(PCMAP) -relief sunken
	place $w.pcme -relwidth .5714 -relheight .0610 -relx .2857 -rely .3341
	helplink $w.pcme trad options pcmap
	button $w.pcmd -text Delete -relief raised -command pcm_delete
	place $w.pcmd -relwidth .1071 -relheight .0610 -relx .98 -rely .3341 \
			-anchor ne
	helplink $w.pcmd trad options pcmdelete
	# Ambfile
	label $w.ambl -text Ambfile:
	place $w.ambl -relx .0714 -rely .4222
	entry $w.ambe -textvariable radvar(AMBFILE) -relief sunken
	place $w.ambe -relwidth .5714 -relheight .0610 -relx .2857 -rely .4222
	helplink $w.ambe trad options ambfile
	button $w.ambd -text Delete -relief raised -command amb_delete
	place $w.ambd -relwidth .1071 -relheight .0610 -relx .98 -rely .4222 \
			-anchor ne
	helplink $w.ambd trad options ambdelete
	# Optfile
	label $w.optl -text Optfile:
	place $w.optl -relx .0714 -rely .5103
	entry $w.opte -textvariable radvar(OPTFILE) -relief sunken
	place $w.opte -relwidth .5714 -relheight .0610 -relx .2857 -rely .5103
	helplink $w.opte trad options optfile
	# Report
	label $w.repl -text Report:
	place $w.repl -relx .0714 -rely .5984
	entry $w.repe -textvariable radvar(REPORT) -relief sunken
	place $w.repe -relwidth .5714 -relheight .0610 -relx .2857 -rely .5984
	helplink $w.repe trad options report
	# Oconv options
	label $w.ocol -text "oconv opts:"
	place $w.ocol -relx .0714 -rely .6863
	entry $w.ocoe -textvariable radvar(oconv) -relief sunken
	place $w.ocoe -relwidth .5714 -relheight .0600 -relx .2857 -rely .6863
	helplink $w.ocoe trad options oconv
	# Mkillum options
	label $w.mkil -text "mkillum opts:"
	place $w.mkil -relx .0714 -rely .7473
	entry $w.mkie -textvariable radvar(mkillum) -relief sunken
	place $w.mkie -relwidth .5714 -relheight .0600 -relx .2857 -rely .7473
	helplink $w.mkie trad options mkillum
	# Mkpmap options
	label $w.mkpl -text "mkpmap opts:"
	place $w.mkpl -relx .0714 -rely .8083
	entry $w.mkpe -textvariable radvar(mkpmap) -relief sunken
	place $w.mkpe -relwidth .5714 -relheight .0600 -relx .2857 -rely .8083
	helplink $w.mkpe trad options mkpmap
	# Render options
	label $w.renl -text "render opts:"
	place $w.renl -relx .0714 -rely .8693
	entry $w.rene -textvariable radvar(render) -relief sunken
	place $w.rene -relwidth .5714 -relheight .0600 -relx .2857 -rely .8693
	helplink $w.rene trad options render
	# Pfilt options
	label $w.pfil -text "pfilt opts:"
	place $w.pfil -relx .0714 -rely .9339
	entry $w.pfie -textvariable radvar(pfilt) -relief sunken
	place $w.pfie -relwidth .5714 -relheight .0600 -relx .2857 -rely .9339
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
