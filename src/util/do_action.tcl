# RCSid: $Id: do_action.tcl,v 2.11 2004/01/01 19:31:46 greg Exp $
#
# Action screen for trad
#

set batch_fmt "Process %d on %s started"
set hostname [exec hostname]

proc make_script {} {		# make run script
	global scname rpview curmess
	set rfn /usr/tmp/rf[pid]
	if {! [save_vars $rfn]} {
		beep
		set curmess "Cannot save variables to temporary file!"
		return
	}
	set radcom {rad -n}
	if {"$rpview" != " ALL"} {
		lappend radcom -v $rpview
	}
	if {"$scname" == {}} {
		set rof /usr/tmp/ro[pid]
	} else {
		set rof $scname
	}
	lappend radcom $rfn > $rof
	set badrun [catch "exec $radcom" curmess]
	if {"$scname" == {}} {
		if {! $badrun} {
			view_txt $rof
		}
		exec rm $rof
	}
	exec rm $rfn
}

proc make_oct args {		# Make octree file ($args is {-t} or {})
	global radvar curmess
	if {"$radvar(scene)" == {} && "$radvar(illum)" == {}} {
		set fi [open "|make $radvar(OCTREE) $args" r]
		while {[gets $fi curmess] >= 0} { update ; after 1000 }
		catch {close $fi} curmess
	}
	if {"$args" == {}} {
		run_rad -v 0
	} else {
		run_rad $args
	}
}

proc run_rad args {		# Run rad command with given arguments
	global curmess
	set rfn /usr/tmp/rf[pid]
	if {! [save_vars $rfn]} {
		beep
		set curmess "Cannot save variables to temporary file!"
		return
	}
	set ot [file mtime $rfn]
	set fi [open "|rad $args $rfn" r]
	while {[gets $fi curmess] >= 0} { update ; after 1000 }
	if {! [catch {close $fi} curmess] && [file mtime $rfn] > $ot} {
		load_vars $rfn {view}
		make_vmenus
	}
	exec rm $rfn
}

proc run_batch {} {		# start rendering in background
	global rpview bfname batch_pid curmess rifname mywin batch_fmt hostname
	if {! [chksave]} {return}
	if [catch {set fo [open $bfname w]} curmess] {
		return
	}
	# Make space for PID to be written later
	puts $fo \
"                                                                             "
	if {$rpview == " ALL"} {
		set radcom "rad"
	} else {
		set radcom "rad -v $rpview"
	}
	lappend radcom "[file tail $rifname]"
	# Put out command
	puts $fo $radcom
	flush $fo
	# Execute in background
	set batch_pid [eval exec $radcom >>& $bfname &]
	# Now, write PID and close file
	seek $fo 0
	puts -nonewline $fo "[format $batch_fmt $batch_pid $hostname]\
				[exec date]"
	close $fo
	set curmess "Batch rendering process $batch_pid started."
	# Correct button stati
	$mywin.rbst configure -state disabled
	$mywin.rbki configure -state normal
	$mywin.rbce configure -state normal
	$mywin.rbvmb configure -state disabled
}

proc kill_batch {} {		# kill batch rendering
	global batch_pid env curmess mywin
	# Kill batch process and children
	switch -glob [exec uname] {
		SunOS -
		ULTRIX {set ps "ps -lg"}
		default {set ps "ps -lu $env(LOGNAME)"}
	}
	set fi [open "|$ps" r]
	gets $fi li
	regexp -indices -nocase {  *PID} $li pidsec
	regexp -indices -nocase {  *PPID} $li ppidsec
	while {[gets $fi li] >= 0} {
		lappend plist "[eval string range {$li} $pidsec]\
				[eval string range {$li} $ppidsec]"
	}
	if [catch {close $fi} curmess] {
		return
	}
	set plist [lsort $plist]
	set procs $batch_pid
	# Two passes required when process id's go into recycling
	for {set i 0} {$i < 2} {incr i} {
		foreach pp $plist {
			if {[lsearch -exact $procs [lindex $pp 1]] >= 0 &&
					($i == 0 || [lsearch -exact $procs \
						[lindex $pp 0]] < 0)} {
				lappend procs [lindex $pp 0]
			}
		}
	}
	if {! [catch {eval exec kill $procs} curmess]} {
		set curmess "Batch process $batch_pid and children killed."
	}
	# change modes
	$mywin.rbki configure -state disabled
	$mywin.rbst configure -state normal
	$mywin.rbvmb configure -state normal
	set batch_pid 0
}

proc make_vmenus {} {		# make/update view menus
	global mywin radvar rvview
	if {"$radvar(view)" == {}} {
		set rvview X
	} else {
		set rvview [lindex [lindex $radvar(view) 0] 0]
	}
	catch {destroy $mywin.rivmb.m}
	catch {destroy $mywin.rbvmb.m}
	if {[llength $radvar(view)] < 2} {
		$mywin.rivmb configure -state disabled
		return
	}
	menu $mywin.rivmb.m
	$mywin.rivmb configure -menu $mywin.rivmb.m -state normal
	foreach v $radvar(view) {
		$mywin.rivmb.m add radiobutton -label [lindex $v 0] \
				-variable rvview -value [lindex $v 0]
	}
	menu $mywin.rbvmb.m
	$mywin.rbvmb configure -menu $mywin.rbvmb.m
	$mywin.rbvmb.m add radiobutton -label ALL \
			-variable rpview -value " ALL"
	$mywin.rbvmb.m add separator
	foreach v $radvar(view) {
		$mywin.rbvmb.m add radiobutton -label [lindex $v 0] \
				-variable rpview -value [lindex $v 0]
	}
}

proc do_action w {		# Action screen
	global rifname rvview rpview radvar bfname batch_pid \
			curmess scname mywin alldone batch_fmt hostname
	if {"$w" == "done"} {
		unset rvview rpview bfname batch_pid mywin
		return
	}
	set bfname [file rootname [file tail $rifname]].err
	frame $w
	set mywin $w
	# Generate octree
	label $w.octl -text {Compile scene only}
	place $w.octl -relx .1429 -rely .0610
	button $w.octb -text oconv -relief raised -command make_oct
	place $w.octb -relwidth .1071 -relheight .0610 -relx .4643 -rely .0610
	helplink $w.octb trad action oconv
	button $w.octt -text Touch -relief raised -command {make_oct -t}
	place $w.octt -relwidth .1071 -relheight .0610 -relx .7143 -rely .0610
	helplink $w.octt trad action touch
	button $w.octf -text Force -relief raised -command \
			{catch {exec rm -f $radvar(OCTREE)} ; make_oct}
	place $w.octf -relwidth .1071 -relheight .0610 -relx .5893 -rely .0610
	helplink $w.octf trad action force
	# Render interactively
	label $w.ril -text {Render interactively}
	place $w.ril -relx .1429 -rely .2439
	button $w.rirv -text rvu -relief raised \
			-command {run_rad -v $rvview -o x11}
	place $w.rirv -relwidth .1071 -relheight .0610 -relx .4643 -rely .2439
	helplink $w.rirv trad action rvu
	label $w.rivl -text View:
	place $w.rivl -relx .6072 -rely .2439
	menubutton $w.rivmb -textvariable rvview -anchor w -relief raised
	place $w.rivmb -relwidth .1071 -relheight .0610 -relx .7143 -rely .2439
	helplink $w.rivmb trad action view
	# Render in background
	label $w.rbl -text {Render in background}
	place $w.rbl -relx .1429 -rely .4268
	button $w.rbst -text Start -relief raised -command run_batch
	place $w.rbst -relwidth .1071 -relheight .0610 -relx .4643 -rely .4268
	helplink $w.rbst trad action start
	label $w.rbvl -text View:
	place $w.rbvl -relx .6072 -rely .4268
	menubutton $w.rbvmb -textvariable rpview -anchor w -relief raised
	place $w.rbvmb -relwidth .1071 -relheight .0610 -relx .7143 -rely .4268
	helplink $w.rbvmb trad action view
	button $w.rbki -text Kill -relief raised -command kill_batch
	place $w.rbki -relwidth .1071 -relheight .0610 -relx .4643 -rely .5488
	helplink $w.rbki trad action kill
	button $w.rbce -text "Check errors" -relief raised \
			-command {view_txt $bfname}
	place $w.rbce -relwidth .1786 -relheight .0610 -relx .6429 -rely .5488
	helplink $w.rbce trad action checkerr
	if [file isfile $bfname] {
		set fi [open $bfname r]
		if {[scan [gets $fi] $batch_fmt batch_pid batch_host] != 2} {
			set batch_host unknown
			set radcom {}
		} else {
			gets $fi radcom
		}
		close $fi
		if [string match "rad -v *" $radcom] {
			set rpview [lindex $radcom 2]
		} else {
			set rpview " ALL"
		}
		if {"$hostname" != "$batch_host"} {
			if $alldone {
				set curmess "Batch rendering finished\
						on $batch_host."
			} else {
				set curmess "Batch rendering on\
						$batch_host -- status UNKNOWN."
				$w.rbst configure -state disabled
				$w.rbvmb configure -state disabled
			}
			set batch_pid 0
			$w.rbki configure -state disabled
		} elseif {[catch {exec /bin/kill -0 $batch_pid}]} {
			if $alldone {
				set curmess "Batch rendering finished."
			} else {
				set curmess "Batch rendering stopped."
			}
			set batch_pid 0
			set rpview " ALL"
			$w.rbki configure -state disabled
		} else {
			set curmess "Batch rendering process $batch_pid running."
			$w.rbst configure -state disabled
			$w.rbvmb configure -state disabled
		}
	} else {
		set curmess "No batch rendering in progress."
		set batch_pid 0
		set rpview " ALL"
		$w.rbki configure -state disabled
		$w.rbce configure -state disabled
	}
	make_vmenus
	# Dry run
	label $w.drl -text {Dry run}
	place $w.drl -relx .1429 -rely .7317
	button $w.drsc -text Script -relief raised -command make_script
	place $w.drsc -relwidth .1071 -relheight .0610 -relx .4643 -rely .7317
	entry $w.drsf -textvariable scname -relief sunken
	place $w.drsf -relwidth .2143 -relheight .0610 -relx .6429 -rely .7317
	helplink "$w.drsc $w.drsf" trad action script
	button $w.drvw -text Edit -relief raised -command {view_txt $scname}
	place $w.drvw -relwidth .1071 -relheight .0610 -relx .4643 -rely .8537
	helplink $w.drvw trad action edit
	button $w.drdel -text Delete -relief raised \
			-command {catch {exec rm $scname < /dev/null} curmess}
	place $w.drdel -relwidth .1071 -relheight .0610 -relx .6429 -rely .8537
	helplink $w.drdel trad action delete
}
