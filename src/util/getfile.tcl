# RCSid: $Id: getfile.tcl,v 2.3 2012/06/03 20:13:51 greg Exp $
#
# Usage: getfile [-win w] [-grab] [-perm] [-glob pattern] [-view proc] [-send proc]
#
# Create a dialog box (in window w if given) to get file name.
# If -grab option is given, then getfile does a local grab on its window.
# Normally, a single file name and return as value.
# If perm switch is given, keep window up for multiple file entries.
# If pattern is given, start with list of all the specified files,
# breaking the string into directory and file match parts.
# If a view procedure is given, a separate "View" button appears
# for allowing the user to preview a selected file.
# If a send procedure is given, a file name is sent upon
# each double-click, and no File entry is displayed.
#

proc getfile args {		# get filename interactively
	# Create necessary globals
	global curdir curfile curpat
	# Set defaults
	set w .fpwin
	set topwin 1
	set dograb 0
	set curdir .
	set curpat *
	set curfile {}
	set transient 1
	# Get options
	while {[llength $args]} {
		switch -glob -- [lindex $args 0] {
			-win* {
				set w [lindex $args 1]
				set topwin 0
				set args [lreplace $args 1 1]
				}
			-grab {
				set dograb 1
				}
			-perm* {
				set transient 0
				}
			-glob {
				set curdir [file dirname [lindex $args 1]]
				set curpat [file tail [lindex $args 1]]
				set args [lreplace $args 1 1]
				}
			-vi* {
				set viewproc [lindex $args 1]
				set args [lreplace $args 1 1]
				}
			-send {
				set sendproc [lindex $args 1]
				set args [lreplace $args 1 1]
				}
			default { error "bad argument to getfile: [lindex $args 0]" }
		}
		set args [lreplace $args 0 0]
	}
	# Save initial directory
	set startdir [pwd]
	# Create widgets
	catch {destroy $w}
	if $topwin {
		toplevel $w -width 400 -height 410
		wm title $w "File Picker"
		wm iconname $w "Files"
		wm minsize $w 400 300
	} else {
		frame $w -width 400 -height 410
		pack $w
	}
	if $dograb { grab $w }
	$w configure -cursor top_left_arrow
	label $w.dt -text Directory:
	place $w.dt -relx .025 -rely .03125
	helplink $w.dt file directory intro
	button $w.ls -text List -command "update_dir $w"
	place $w.ls -relx .025 -rely .125 -relwidth .15 -relheight .0625
	helplink $w.ls file directory list
	entry $w.de -textvariable curdir -relief sunken
	place $w.de -relx .25 -rely .03125 -relwidth .6875 -relheight .0625
	focus $w.de
	helplink $w.de file directory intro
	bind $w.de <Return> "update_dir $w"
	entry $w.pw -relief sunken -textvariable curpat
	place $w.pw -relx .25 -rely .125 -relwidth .34375 -relheight .0625
	bind $w.pw <Return> "update_dir $w"
	helplink $w.pw file directory match
	frame $w.fm
	scrollbar $w.fm.sb -relief sunken -command "$w.fm.fl yview"
	listbox $w.fm.fl -relief sunken -yscroll "$w.fm.sb set" \
			-selectmode browse
	pack $w.fm.sb -side right -fill y
	pack $w.fm.fl -side left -expand yes -fill both
	place $w.fm -relx .25 -rely .21875 -relwidth .6875 -relheight .625
	helplink $w.fm.fl file file intro
	if [info exists sendproc] {
		bind $w.fm.fl <Double-Button-1> "set_curfile $w $sendproc"
		bind $w.fm.fl <Return> "set_curfile $w $sendproc"
		if {$transient} {
			bind $w.fm.fl <Double-Button-1> \
					"+if {\$curfile > {}} {destroy $w}"
			bind $w.fm.fl <Return> \
					"+if {\$curfile > {}} {destroy $w}"
		}
	} else {
		bind $w.fm.fl <Double-Button-1> "set_curfile $w"
		bind $w.fm.fl <Return> "set_curfile $w"
		label $w.fl -text File:
		place $w.fl -relx .10625 -rely .875
		entry $w.fn -relief sunken -textvariable curfile
		place $w.fn -relx .25 -rely .875 \
				-relwidth .6875 -relheight .0625
		focus $w.fn
		helplink $w.fn file file entry
		if $transient {
			bind $w.fn <Return> "destroy $w"
		} else {
			bind $w.fn <Return> {}
		}
	}
	if {$transient != [info exists sendproc]} {
		button $w.ok -text OK -command "destroy $w"
		place $w.ok -relx .025 -rely .28125 \
				-relwidth .15 -relheight .0625
	}
	if {$transient || [info exists sendproc]} {
		button $w.cancel -text Cancel \
				-command "destroy $w; unset curdir"
		place $w.cancel -relx .025 -rely .375 \
				-relwidth .15 -relheight .0625
	}
	if [info exists viewproc] {
		button $w.vi -text View -state disabled \
		-command "$viewproc \[$w.fm.fl get \[$w.fm.fl curselection\]\]"
		place $w.vi -relx .025 -rely .46875 \
				-relwidth .15 -relheight .0625
		bind $w.fm.fl <ButtonRelease-1> "+chk_select $w"
		helplink $w.vi file file view
	}
	# Load current directory
	update_dir $w
	if $transient {
		tkwait window $w
		cd $startdir
		if [info exists sendproc] {
			return [info exists curdir]
		}
		if [info exists curdir] {
			return $curdir/$curfile
		}
		return
	} elseif {[info exists sendproc]} {
		tkwait window $w
		cd $startdir
		if [info exists curdir] {
			return 1
		}
		return 0
	}
	# Else return the dialog window
	return $w
}


proc update_dir w {			# Update working directory
	global curdir curpat
	if {"$curpat" == {}} {
		set curpat *
	}
	if {"$curdir" == {}} {set curdir {~}}
	if [catch {cd $curdir} curdir] {
		set oldcol [lindex [$w.de config -foreground] 4]
		$w.de config -foreground Red
		update idletasks
		after 1500
		$w.de config -foreground $oldcol
	}
	set curdir [pwd]
	set offset [expr [string length $curdir] - 32]
	if {$offset > 0} {
		$w.de xview [expr $offset + \
			[string first / [string range $curdir $offset end]] \
				+ 1]
	}
	$w.fm.fl delete 0 end
	set ls ../
	foreach f [glob -nocomplain *] {
		if [file isdirectory $f] {
			lappend ls $f/
		}
	}
	foreach f [eval glob -nocomplain [split $curpat]] {
		if [file isfile $f] {
			lappend ls $f
		}
	}
	eval $w.fm.fl insert end [lsort $ls]
}


proc set_curfile {w {sp {}}} {		# change current file selection
	global curdir curfile
	set f [$w.fm.fl get [$w.fm.fl curselection]]
	if [string match */ $f] {
		set curdir [string range $f 0 [expr [string length $f]-2]]
		update_dir $w
	} else {
		set curfile $f
		if {"$sp" != {}} {
			$sp $curdir/$curfile
		}
	}
}


proc chk_select w {			# check if current selection is file
	set s [$w.fm.fl curselection]
	if {"$s" != {} && [file isfile [$w.fm.fl get $s]]} {
		$w.vi configure -state normal
	} else {
		$w.vi configure -state disabled
	}
}
