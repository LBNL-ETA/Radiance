# RCSid: $Id: do_file.tcl,v 2.18 2005/02/16 05:40:12 greg Exp $
#
# Choose the Rad Input File to work on.
#

proc preen {} {			# clean up radvar
	global radvar rifname
	foreach n {objects scene materials illum mkillum render oconv pfilt
			RAWFILE ZFILE AMBFILE OPTFILE EXPOSURE ZONE REPORT} {
		if {! [info exists radvar($n)]} {
			set radvar($n) {}
		}
	}
	if [info exists radvar(view)] {
		set oldval $radvar(view)
		set radvar(view) {}
		set n 1
		foreach v $oldval {
			if {"[string index $v 0]" == "-"} {
				lappend radvar(view) "u$n $v"
			} elseif {[lsearch -glob $radvar(view) \
					"[lindex $v 0] *"] >= 0} {
				continue
			} else {
				lappend radvar(view) $v
			}
			incr n
		}
	} else {
		set radvar(view) {}
	}
	if {! [info exists radvar(UP)]} {
		set radvar(UP) Z
	} elseif {[string match +* $radvar(UP)]} {
		set radvar(UP) [string toupper [string range $radvar(UP) 1 end]]
	} else {
		set radvar(UP) [string toupper $radvar(UP)]
	}
	set rifroot [file rootname [file tail $rifname]]
	if {! [info exists radvar(OCTREE)]} {
		set radvar(OCTREE) $rifroot.oct
	}
	if {! [info exists radvar(RESOLUTION)]} {
		set radvar(RESOLUTION) 512
	}
	if {! [info exists radvar(EYESEP)]} {
		set radvar(EYESEP) 1
	}
	if [info exists radvar(QUALITY)] {
		cardval radvar(QUALITY) {High Medium Low}
	} else {
		set radvar(QUALITY) Low
	}
	if {! [info exists radvar(PICTURE)]} {
		set radvar(PICTURE) $rifroot
	}
	if {! [info exists radvar(INDIRECT)]} {
		set radvar(INDIRECT) 0
	}
	if [info exists radvar(DETAIL)] {
		cardval radvar(DETAIL) {High Medium Low}
	} else {
		set radvar(DETAIL) Medium
	}
	if [info exists radvar(PENUMBRAS)] {
		cardval radvar(PENUMBRAS) {True False}
	} else {
		set radvar(PENUMBRAS) False
	}
	if [info exists radvar(VARIABILITY)] {
		cardval radvar(VARIABILITY) {High Medium Low}
	} else {
		set radvar(VARIABILITY) Low
	}
}

proc setradvar stmt {		# assign a rad variable
	global radvar
	regexp {^([a-zA-Z][a-zA-Z0-9]*) *=[ 	]*(.*)$} $stmt dummy vnam vval
	switch -glob $vnam {
		obj* { eval lappend radvar(objects) $vval }
		sce* { eval lappend radvar(scene) $vval }
		mat* { eval lappend radvar(materials) $vval }
		ill* { eval lappend radvar(illum) $vval }
		mki* { eval lappend radvar(mkillum) $vval }
		ren* { eval lappend radvar(render) $vval }
		oco* { eval lappend radvar(oconv) $vval }
		pf* { eval lappend radvar(pfilt) $vval }
		vi* { lappend radvar(view) $vval }
		ZO* { set radvar(ZONE) $vval }
		QUA* { set radvar(QUALITY) $vval }
		OCT* { set radvar(OCTREE) $vval }
		PIC* { set radvar(PICTURE) $vval }
		AMB* { set radvar(AMBFILE) $vval }
		OPT* { set radvar(OPTFILE) $vval }
		EXP* { set radvar(EXPOSURE) $vval }
		EYE* { set radvar(EYESEP) $vval }
		RES* { set radvar(RESOLUTION) $vval }
		UP { set radvar(UP) $vval }
		IND* { set radvar(INDIRECT) $vval }
		DET* { set radvar(DETAIL) $vval }
		PEN* { set radvar(PENUMBRAS) $vval }
		VAR* { set radvar(VARIABILITY) $vval }
		REP* { set radvar(REPORT) $vval }
		RAW* { set radvar(RAWFILE) $vval }
		ZF* {set radvar(ZFILE) $vval }
	}
		
}

proc putradvar {fi vn} {	# print out a rad variable
	global radvar
	if {! [info exists radvar($vn)] || $radvar($vn) == {}} {return}
	if [regexp {^[A-Z]} $vn] {
		puts $fi "$vn= $radvar($vn)"
		return
	}
	if {"$vn" == "view"} {
		foreach v $radvar(view) {
			puts $fi "view= $v"
		}
		return
	}
	if {[lsearch -exact {ZONE QUALITY OCTREE PICTURE AMBFILE OPTFILE
			EXPOSURE RESOLUTION UP INDIRECT DETAIL PENUMBRAS
			EYESEP RAWFILE ZFILE VARIABILITY REPORT} $vn] >= 0} {
		puts $fi "$vn= $radvar($vn)"
		return
	}
	puts -nonewline $fi "$vn="
	set vnl [expr [string length $vn] + 1]
	set pos $vnl
	for {set i 0} {$i < [llength $radvar($vn)]} {incr i} {
		set len [expr [string length [lindex $radvar($vn) $i]] + 1]
		if {$pos > $vnl && $pos + $len > 70} {
			puts -nonewline $fi "\n$vn="
			set pos $vnl
		}
		puts -nonewline $fi " [lindex $radvar($vn) $i]"
		incr pos $len
	}
	puts $fi {}
}

proc load_vars {f {vl all}} {	# load RIF variables
	global curmess radvar alldone
	if {"$f" == ""} {return 0}
	if {! [file isfile $f]} {
		beep
		set curmess "$f: no such file."
		return 0
	}
	if {"$vl" == "all" && ! [chksave]} {return 0}
	set curmess {Please wait...}
	update
	if [catch {exec rad -n -w -e $f >& /tmp/ro[pid]}] {
		set curmess [exec cat /tmp/ro[pid]]
		exec rm -f /tmp/ro[pid]
		return 0
	}
	set fi [open /tmp/ro[pid] r]
	if {"$vl" == "all"} {
		catch {unset radvar}
		while {[gets $fi curli] != -1} {
			if [regexp {^[a-zA-Z][a-zA-Z0-9]* *=} $curli] {
				setradvar $curli
			} else {
				break
			}
		}
		set curmess {Project loaded.}
	} else {
		foreach n $vl {
			if [regexp {[a-z][a-z]*} $n] {
				set radvar($n) {}
			} else {
				catch {unset radvar($n)}
			}
		}
		while {[gets $fi curli] != -1} {
			if [regexp {^[a-zA-Z][a-zA-Z0-9]* *=} $curli] {
				regexp {^[a-zA-Z][a-zA-Z0-9]*} $curli thisv
				if {[lsearch -exact $vl $thisv] >= 0} {
					setradvar $curli
				}
			} else {
				break
			}
		}
		set curmess {Variables loaded.}
	}
	set alldone [eof $fi]
	close $fi
	exec rm -f /tmp/ro[pid]
	preen
	return 1
}

proc save_vars {f {vl all}} {	# save RIF variables
	global curmess radvar
	if {"$f" == {} || ! [info exists radvar]} {return 0}
	if {"$vl" == "all"} {
		if [catch {set fi [open $f w]} curmess] {
			beep
			return 0
		}
		puts $fi "# Rad Input File created by trad [exec date]"
		foreach n [lsort [array names radvar]] {
			putradvar $fi $n
		}
	} else {
		if [catch {set fi [open $f a]} curmess] {
			beep
			return 0
		}
		foreach n [array names radvar] {
			if {[lsearch -exact $vl $n] >= 0} {
				putradvar $fi $n
			}
		}
	}
	close $fi
	set curmess {File saved.}
	return 1
}

proc newload f {		# load a RIF
	global rifname readonly
	if [load_vars $f] {
		set rifname [pwd]/$f
		set readonly [expr ! [file writable $f]]
		gotfile 1
		return 1
	}
	return 0
}

proc newsave f {		# save a RIF
	global rifname readonly curmess
	if {"[pwd]/$f" == "$rifname"} {
		if $readonly {
			beep
			set curmess {Mode set to read only.}
			return 0
		}
	} elseif {[file exists $f]} {
		set ftyp [file type $f]
		if { $ftyp != "file" } {
			beep
			set curmess "Selected file $f is a $ftyp."
			return 0
		}
		if [tk_dialog .dlg {Verification} \
				"Overwrite existing file $f?" \
				questhead 1 {Go Ahead} {Cancel}] {
			return 0
		}
	}
	if {[file isfile $f] && ! [file writable $f] && 
			[catch {exec chmod u+w $f} curmess]} {
		beep
		return 0
	}
	if [save_vars $f] {
		set rifname [pwd]/$f
		set readonly 0
		gotfile 1
		return 1
	}
	return 0
}

proc newnew f {			# create a new RIF
	global rifname readonly curmess radvar
	if [file exists $f] {
		set ftyp [file type $f]
		if { $ftyp != "file" } {
			beep
			set curmess "Selected file $f is a $ftyp."
			return 0
		}
		if [tk_dialog .dlg {Verification} \
				"File $f exists -- disregard it?" \
				questhead 1 {Yes} {Cancel}] {
			return 0
		}
	}
	if {! [chksave]} {return 0}
	set rifname [pwd]/$f
	set readonly 0
	catch {unset radvar}
	preen
	gotfile 1
	set curmess {Go to SCENE and enter OCTREE and/or scene file(s).}
	return 1
}

proc do_file w {
	global rifname readonly myglob curfile curpat
	if {"$w" == "done"} {
		cd [file dirname $rifname]
		set myglob(rif) $curpat
		return
	}
	frame $w
	frame $w.left
	pack $w.left -side left
	button $w.left.load -text LOAD -width 5 \
			-relief raised -command {newload $curfile}
	button $w.left.save -text SAVE -width 5 -relief raised \
			-command "newsave \$curfile; update_dir $w.right"
	button $w.left.new -text NEW -width 5 \
			-relief raised -command {newnew $curfile}
	pack $w.left.load $w.left.save $w.left.new -side top -pady 15 -padx 20
	checkbutton $w.left.ro -text "Read Only" -variable readonly \
			-onvalue 1 -offvalue 0 -relief flat
	pack $w.left.ro -side top -pady 15 -padx 20
	helplink $w.left.load trad file load
	helplink $w.left.save trad file save
	helplink $w.left.new trad file new
	helplink $w.left.ro trad file readonly
	getfile -view view_txt -perm \
			-win $w.right -glob [file dirname $rifname]/$myglob(rif)
	set curfile [file tail $rifname]
}
