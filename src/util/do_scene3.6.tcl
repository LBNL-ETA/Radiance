# RCSid: $Id: do_scene3.6.tcl,v 2.10 2003/02/22 02:07:30 greg Exp $
#
# Track octree and scene files
#

proc newfent f {		# add file to our list
	global rifname radvar mybox myvar
	set rd [file dirname $rifname]
	if {[string first $rd $f] == 0} {
		set f [string range $f [expr [string length $rd] + 1] end]
	}
	if {[lsearch -exact $radvar($myvar) $f] == -1} {
		lappend radvar($myvar) $f
		$mybox($myvar) insert end $f
		$mybox($myvar) yview end
	}
}

proc lbgetf nm {		# get list box files
	global myglob radvar mybox myvar curpat curmess
	set myvar $nm
	set oldnum [llength $radvar($nm)]
	if [getfile -grab -perm -glob $myglob($nm) -view view_txt -send newfent] {
		set curmess "Added [expr [llength $radvar($nm)] - $oldnum] entries."
	} elseif {[llength $radvar($nm)] > $oldnum} {
		set radvar($nm) [lreplace $radvar($nm) $oldnum end]
		set curmess "Cancelled."
	}
	$mybox($nm) delete 0 end
	eval $mybox($nm) insert end $radvar($nm)
	set myglob($nm) $curpat
	unset myvar
}

proc oct_delete {} {		# delete octree file
	global radvar curmess
	if {"$radvar(OCTREE)" == {} || ! [file isfile $radvar(OCTREE)]} {
		set curmess {No octree file.}
		return
	}
	if [tk_dialog .dlg {Verification} \
			"Really delete octree file $radvar(OCTREE)?" \
			questhead 0 {Delete} {Cancel}] {
		return 0
	}
	if [catch {exec rm $radvar(OCTREE) < /dev/null} curmess] {return 0}
	set curmess {Octree file deleted.}
	return 1
}

proc getdepend {} {		# get object dependencies
	global radvar curmess mybox
	set curmess "Please wait..."
	update
	foreach newf [eval exec raddepend $radvar(illum) $radvar(scene)] {
		if {[lsearch -exact $radvar(objects) $newf] < 0} {
			lappend radvar(objects) $newf
		}
	}
	$mybox(objects) delete 0 end
	eval $mybox(objects) insert end $radvar(objects)
	set curmess "Done."
}

proc vwselfil {} {		# View selected file entries
	global mybox
	foreach n {materials illum scene objects} {
		foreach sl [$mybox($n) curselection] {
			lappend files [$mybox($n) get $sl]
		}
	}
	view_txt $files
}

proc movselfil {n y} {		# move selected files to new position
	global radvar mybox curmess
	set dl [delselfil]
	if {"$dl" == {}} {		# get it from another window
		set dl [selection get]
		set curmess "Pasted [llength $dl] entries."
	} else {
		set curmess "Moved [llength $dl] entries."
	}
	if [llength $dl] {
		# The following should return "end" if past end, but doesn't!
		set i [$mybox($n) nearest $y]
		# So, we hack rather badly to discover the truth...
		if {$i == [$mybox($n) size] - 1 && $y > 12 &&
				[$mybox($n) nearest [expr $y - 12]] == $i} {
			set i end
		}
		eval $mybox($n) insert $i $dl
		if {"$i" == "end"} {
			eval lappend radvar($n) $dl
		} else {
			set radvar($n) [eval linsert {$radvar($n)} $i $dl]
		}
	}
}

proc delselfil {} {		# Delete selected file entries
	global radvar mybox curmess
	set dl {}
	foreach n {materials illum scene objects} {
		foreach sl [lsort -integer -decreasing [$mybox($n) curselection]] {
			set i [lsearch -exact $radvar($n) [$mybox($n) get $sl]]
			set dl "[lindex $radvar($n) $i] $dl"
			set radvar($n) [lreplace $radvar($n) $i $i]
			$mybox($n) delete $sl
		}
	}
	set curmess "Discarded [llength $dl] entries."
	return $dl
}

proc copyscene rf {		# Copy scene data from specified RIF
	global mybox radvar
	load_vars [file tail $rf] {OCTREE materials illum scene objects}
	foreach n {materials illum scene objects} {
		 $mybox($n) delete 0 end
		eval $mybox($n) insert end $radvar($n)
	}
}

proc do_scene w {		# Create scene screen
	global radvar mybox rifname
	if {"$w" == "done"} {
		unset mybox
		return
	}
	set lbfont -*-courier-medium-r-normal--12-*-*-*-*-*-iso8859-1
	frame $w
	# Octree entry
	label $w.octl -text Octree
	place $w.octl -relwidth .1071 -relheight .0610 -relx .0714 -rely .0610
	entry $w.oct -textvariable radvar(OCTREE) -relief sunken
	place $w.oct -relwidth .5714 -relheight .0610 -relx .2143 -rely .0610
	helplink $w.oct trad scene octree
	button $w.odel -text Delete -relief raised -command oct_delete
	place $w.odel -relwidth .1071 -relheight .0610 -relx .8000 -rely .0610
	helplink $w.odel trad scene octdelete
	# Materials listbox
	button $w.matb -text Materials -relief raised \
			-command "lbgetf materials"
	set mybox(materials) $w.mat.lb
	place $w.matb -relwidth .1200 -relheight .0610 -relx .0714 -rely .1463
	frame $w.mat
	scrollbar $w.mat.sb -relief sunken -command "$w.mat.lb yview"
	listbox $w.mat.lb -relief sunken -yscroll "$w.mat.sb set" -font $lbfont
	bind $w.mat.lb <Button-2> "movselfil materials %y"
	pack $w.mat.sb -side right -fill y
	pack $w.mat.lb -side left -expand yes -fill both
	place $w.mat -relwidth .5714 -relheight .0976 -relx .2143 -rely .1463
	eval $w.mat.lb insert end $radvar(materials)
	helplink "$w.mat.lb $w.matb" trad scene materials
	# Illum listbox
	button $w.illb -text Illum -relief raised \
			-command "lbgetf illum"
	set mybox(illum) $w.ill.lb
	place $w.illb -relwidth .1200 -relheight .0610 -relx .0714 -rely .2683
	frame $w.ill
	scrollbar $w.ill.sb -relief sunken -command "$w.ill.lb yview"
	listbox $w.ill.lb -relief sunken -yscroll "$w.ill.sb set" -font $lbfont
	bind $w.ill.lb <Button-2> "movselfil illum %y"
	pack $w.ill.sb -side right -fill y
	pack $w.ill.lb -side left -expand yes -fill both
	place $w.ill -relwidth .5714 -relheight .0976 -relx .2143 -rely .2683
	eval $w.ill.lb insert end $radvar(illum)
	helplink "$w.ill.lb $w.illb" trad scene illum
	# Scene listbox
	button $w.sceb -text Scene -relief raised \
			-command "lbgetf scene"
	set mybox(scene) $w.sce.lb
	place $w.sceb -relwidth .1200 -relheight .0610 -relx .0714 -rely .3902
	frame $w.sce
	scrollbar $w.sce.sb -relief sunken -command "$w.sce.lb yview"
	listbox $w.sce.lb -relief sunken -yscroll "$w.sce.sb set" -font $lbfont
	bind $w.sce.lb <Button-2> "movselfil scene %y"
	pack $w.sce.sb -side right -fill y
	pack $w.sce.lb -side left -expand yes -fill both
	place $w.sce -relwidth .5714 -relheight .2683 -relx .2143 -rely .3902
	eval $w.sce.lb insert end $radvar(scene)
	helplink "$w.sce.lb $w.sceb" trad scene scene
	# Objects listbox
	button $w.objb -text Objects -relief raised \
			-command "lbgetf objects"
	set mybox(objects) $w.obj.lb
	place $w.objb -relwidth .1200 -relheight .0610 -relx .0714 -rely .6829
	frame $w.obj
	scrollbar $w.obj.sb -relief sunken -command "$w.obj.lb yview"
	listbox $w.obj.lb -relief sunken -yscroll "$w.obj.sb set" -font $lbfont
	bind $w.obj.lb <Button-2> "movselfil objects %y"
	pack $w.obj.sb -side right -fill y
	pack $w.obj.lb -side left -expand yes -fill both
	place $w.obj -relwidth .5714 -relheight .2683 -relx .2143 -rely .6829
	eval $w.obj.lb insert end $radvar(objects)
	button $w.autob -text Auto -relief raised -command getdepend
	place $w.autob -relwidth .1200 -relheight .0610 -relx .0714 -rely .7927
	helplink "$w.obj.lb $w.objb $w.autob" trad scene objects
	# View button
	button $w.vwb -text Edit -relief raised -command vwselfil
	place $w.vwb -relwidth .1071 -relheight .0610 -relx .8214 -rely .4000
	helplink $w.vwb trad scene edit
	# Delete button
	button $w.del -text Discard -relief raised -command delselfil
	place $w.del -relwidth .1071 -relheight .0610 -relx .8214 -rely .5000
	helplink $w.del trad scene discard
	# Revert and Copy buttons
	button $w.revert -text Revert -relief raised \
			-command "copyscene $rifname"
	place $w.revert -relwidth .1071 -relheight .0610 -relx .98 -rely .98 \
			-anchor se
	helplink $w.revert trad scene revert
	button $w.copy -text Copy -relief raised -command {getfile -grab \
			-send copyscene -view view_txt -glob $rif_glob}
	place $w.copy -relwidth .1071 -relheight .0610 -relx .98 -rely .90 \
			-anchor se
	helplink $w.copy trad scene copy
}
