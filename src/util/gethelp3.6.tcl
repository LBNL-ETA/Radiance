# RCSid: $Id: gethelp3.6.tcl,v 2.12 2003/02/22 02:07:30 greg Exp $
#
# Get help from a file, formatted like so:
#
# .category1.topic1
#
# A bunch of unformatted text...
#
# .category2.topic2
#
# etc.
#

# Calling program should define global help file directory as $helplib

proc helplink {win file cat top} {		# Create link for window(s)
	foreach w $win {
		bind $w <Control-Button-1> "gethelp $file $cat $top"
	}
}

set curhelp(file) {}

proc gethelp {helpfile category topic} {	# Open help window
	global curhelp helpfontwidth
	if {! [winfo exists .helpwin]} {		# Set up window
		toplevel .helpwin -cursor top_left_arrow
		wm minsize .helpwin 500 400
		wm iconbitmap .helpwin question
		frame .helpwin.but -width 150 -height 400
		pack .helpwin.but -side right -fill none -expand no
		label .helpwin.but.lab -textvariable curhelp(title)
		place .helpwin.but.lab -relx .1667 -rely 0
		helplink .helpwin.but.lab help help intro
		menubutton .helpwin.but.catb -relief raised -height 1 \
				-text Category -menu .helpwin.but.catb.m
		menu .helpwin.but.catb.m
		place .helpwin.but.catb -relx .1667 -rely .1100 \
				-relwidth .5000 -relheight .0600
		helplink .helpwin.but.catb help help category
		menubutton .helpwin.but.topb -relief raised -height 1 \
				-text Topic -menu .helpwin.but.topb.m
		menu .helpwin.but.topb.m
		place .helpwin.but.topb -relx .1667 -rely .2200 \
				-relwidth .5000 -relheight .0600
		helplink .helpwin.but.topb help help topic
		button .helpwin.but.srchb -text Grep -relief raised \
				-height 1 -command {helpsearch $curhelp(search)}
		place .helpwin.but.srchb -relx .1667 -rely .3500 \
				-relwidth .5000 -relheight .0600
		helplink .helpwin.but.srchb help navigate search
		entry .helpwin.but.srche -relief sunken -insertofftime 0 \
				-textvariable curhelp(search)
		place .helpwin.but.srche -relx .1667 -rely .4200 \
				-relwidth .6667 -relheight .0600
		bind .helpwin.but.srche <Return> {set curhelp(msg) {}
							helpupdate}
		helplink .helpwin.but.srche help navigate search
		label .helpwin.but.msg -textvariable curhelp(msg)
		place .helpwin.but.msg -relx .1667 -rely .4800
		button .helpwin.but.next -text Next -relief raised -height 1 \
				-command {eval helphist new $curhelp(next)}
		place .helpwin.but.next -relx .1667 -rely .7500 \
				-relwidth .5 -relheight .0600
		helplink .helpwin.but.next help navigate next
		button .helpwin.but.back -text Back -relief raised -height 1 \
				-command "helphist back"
		place .helpwin.but.back -relx .1667 -rely .6500 \
				-relwidth .5 -relheight .0600
		helplink .helpwin.but.back help navigate back
		button .helpwin.but.forw -text Forward -relief raised -height 1\
				-command "helphist forward"
		place .helpwin.but.forw -relx .1667 -rely .5500 \
				-relwidth .5 -relheight .0600
		helplink .helpwin.but.forw help navigate forward
		button .helpwin.but.done -text Done -relief raised -height 1 \
				-command "destroy .helpwin ; helpopen {}"
		place .helpwin.but.done -relx .1667 -rely .9000 \
				-relwidth .5 -relheight .0600
		helplink .helpwin.but.done help help done
		text .helpwin.txt -wrap word -width 48 -height 20 -bd 2 \
			-yscrollcommand ".helpwin.sb set" -relief raised \
			-font -*-courier-medium-r-normal--14-*-*-*-*-*-iso8859-1
		.helpwin.txt tag configure highlight \
			-font -*-courier-bold-r-normal--14-*-*-*-*-*-iso8859-1
		set helpfontwidth 9
		scrollbar .helpwin.sb -relief flat \
				-command ".helpwin.txt yview"
		pack .helpwin.sb -side right -fill y
		pack .helpwin.txt -expand yes -fill both
		helplink .helpwin.sb help helpwin scroll
		helplink .helpwin.txt help helpwin intro
		tkwait visibility .helpwin
	} elseif {! [winfo ismapped .helpwin]} {	# map window
		wm deiconify .helpwin
	} else {					# raise window
		raise .helpwin
	}
	focus .helpwin.but.srche
	helpopen $helpfile
	helphist new $category $topic
}

proc helpopen fname {			# open the named help file
	global curhelp helpindex helplib
	if {"$fname" == "$curhelp(file)"} {return}
	if {"$curhelp(file)" != {}} {
		close $curhelp(fid)
		unset curhelp helpindex
	}
	if {[set curhelp(file) $fname] == {}} {return}
	# Complete file name as required
	if {! [info exists helplib]} {
		set helplib /usr/local/lib/tk/help
	}
	if {"[file rootname $fname]" == "$fname"} {
		append fname .hlp
	}
	if {! [string match {[~/.]*} $fname]} {
		set fname $helplib/$fname
	}
	wm title .helpwin $fname
	set curhelp(title) "[string toupper\
			[file rootname [file tail $fname]]]  HELP"
	set ifile [file rootname $fname].ndx
	wm iconname .helpwin [string tolower $curhelp(title)]
	.helpwin.txt configure -state normal
	.helpwin.txt delete 1.0 end
	.helpwin.txt insert end "Loading $fname..."
	update
	set curhelp(fid) [open $fname r]
	if {! [file isfile $ifile] ||
			[file mtime $fname] > [file mtime $ifile]} {
		set helpindex(catlist) {}
		while {[gets $curhelp(fid) li] >= 0} {
			if [regexp -nocase {^\.([A-Z][A-Z0-9]*)\.([A-Z][A-Z0-9]*)$} \
					$li dummy cat top] {
				lappend helpindex([string toupper $cat]) $top
				set helpindex([string toupper $cat,$top]) \
						[tell $curhelp(fid)]
				if {[lsearch -exact $helpindex(catlist) $cat] < 0} {
					lappend helpindex(catlist) $cat
				}
			}
		}
		if {! [catch {set fi [open $ifile w]}]} {
			puts $fi "# This is an automatically created index\
					file -- DO NOT EDIT!"
			writevars $fi helpindex
			close $fi
			catch {exec chmod 666 $ifile}
		}
	} else {
		source $ifile
	}
	.helpwin.but.catb.m delete 0 last
	foreach cat $helpindex(catlist) {
		.helpwin.but.catb.m add command -label $cat \
				-command "helphist new $cat intro"
	}
	set curhelp(category) None
	.helpwin.but.topb.m delete 0 last
	set curhelp(topic) None
	set curhelp(next) {}
	helphist clear
}

proc helpgoto {cat top} {		# find selected category and topic
	global curhelp helpindex
	# Capitalize "just in case"
	set cat [string toupper $cat]
	set top [string toupper $top]
	# Change topic menu if category is changed
	set curhelp(topic) $top
	if {"$cat" != "$curhelp(category)"} {
		set curhelp(category) $cat
		.helpwin.but.topb.m delete 0 last
		foreach top $helpindex($cat) {
			.helpwin.but.topb.m add command -label $top \
				-command "helphist new \$curhelp(category) $top"
		}
	}
	helpupdate
}

proc helpupdate {} {			# update help text window
	global curhelp helpindex helpfontwidth
	# Print category and topic
	set linelen [expr "[winfo width .helpwin.txt] / $helpfontwidth - 1"]
	.helpwin.txt configure -state normal
	.helpwin.txt delete 1.0 end
	set titlen [expr [string length $curhelp(category)] + \
			[string length $curhelp(topic)] + 1]
	for {set i [expr ($linelen - $titlen)/2]} {[incr i -1] >= 0} {} {
		.helpwin.txt insert end { }
	}
	.helpwin.txt insert end "$curhelp(category) $curhelp(topic)\n"
	.helpwin.txt tag add highlight "1.0 lineend - $titlen c" "1.0 lineend"
	# Search for it in file
	if {"$curhelp(category) $curhelp(topic)" !=
			"[string toupper $curhelp(next)]"} {
		if [info exists helpindex($curhelp(category),$curhelp(topic))] {
			seek $curhelp(fid) $helpindex($curhelp(category),$curhelp(topic))
		} else {
			.helpwin.txt insert end "\nNo such help topic."
			.helpwin.txt configure -state disabled
			return ".$curhelp(category).$curhelp(topic) not found\
					in $curhelp(file)"
		}
	}
	# Load help text into our window
	set linepos 0
	while {[set ll [gets $curhelp(fid) li]] >= 0 && ! [regexp -nocase \
			{^\.([A-Z][A-Z0-9]*)\.([A-Z][A-Z0-9]*)$} \
			$li dummy cat top]} {
		if {$ll == 0} {				# paragraph
			if $linepos {
				.helpwin.txt insert end "\n\n"
				set linepos 0
			} else {
				.helpwin.txt insert end "\n"
			}
		} else {				# line
			.helpwin.txt insert end $li
			incr linepos $ll
			.helpwin.txt insert end { }
			incr linepos
			# Highlight search string match
			if {"$curhelp(search)" != {}} {
				if [regexp -nocase -indices \
						$curhelp(search) $li mi] {
					.helpwin.txt tag add highlight\
							"end - 1 c - $ll c\
							+ [lindex $mi 0] c"\
							"end - $ll c\
							+ [lindex $mi 1] c"
				}
			}
			# Add extra space at the end of a sentence
			if [regexp {[.?!:][)"']?$} $li] {
				.helpwin.txt insert end { }
				incr linepos
			}
		}
	}
	# Highlight next category and topic
	if {$ll > 0} {
		.helpwin.txt insert end "Next:  $cat $top"
		.helpwin.txt tag add highlight "end linestart" end
		set curhelp(next) "$cat $top"
		.helpwin.but.next configure -state normal
	} else {
		set curhelp(next) {}
		.helpwin.but.next configure -state disabled
	}
	.helpwin.txt configure -state disabled
	return {}
}

proc helphist {op args} {		# access help history list
	global curhelp
	switch -exact $op {
		clear {
			set curhelp(histlst) {}
			set curhelp(histpos) 0
			set gotoit 0
		}
		new {
			set curhelp(histpos) \
					[expr [llength $curhelp(histlst)] - 1]
			if {"[string toupper [lindex $curhelp(histlst) \
					$curhelp(histpos)]]"
					!= "[string toupper $args]"} {
				lappend curhelp(histlst) $args
				incr curhelp(histpos)
			}
			set gotoit 1
		}
		append {
			lappend curhelp(histlst) $args
			set gotoit 0
		}
		back {
			incr curhelp(histpos) -1
			set gotoit 1
		}
		forward {
			incr curhelp(histpos) 1
			set gotoit 1
		}
	}
	if {$curhelp(histpos) + 1 >= [llength $curhelp(histlst)]} {
		.helpwin.but.forw configure -state disabled
	} else {
		.helpwin.but.forw configure -state normal
	}
	if {$curhelp(histpos) <= 0} {
		.helpwin.but.back configure -state disabled
	} else {
		.helpwin.but.back configure -state normal
	}
	if $gotoit {
		eval helpgoto [lindex $curhelp(histlst) $curhelp(histpos)]
	}
}

proc helpsearch word {		# search for occurances of the given word
	global curhelp helpindex
	if {"$curhelp(search)" == {}} {
		set curhelp(msg) "No pattern."
		return 0
	}
	set nmatches 0
	set cat [lindex $helpindex(catlist) 0]
	set top [lindex $helpindex([string toupper $cat]) 0]
	set startpos [tell $curhelp(fid)]
	seek $curhelp(fid) $helpindex([string toupper $cat,$top])
	set foundmatch 0
	while {[gets $curhelp(fid) li] >= 0} {
		if [regexp -nocase {^\.([A-Z][A-Z0-9]*)\.([A-Z][A-Z0-9]*)$} \
				$li dummy cat top] {
			set foundmatch 0
		} elseif {! $foundmatch && \
				[regexp -nocase $curhelp(search) $li]} {
			set foundmatch 1
			if $nmatches {
				helphist append $cat $top
			} else {
				helphist new $cat $top
			}
			incr nmatches
		}
	}
	if $nmatches {
		set curhelp(msg) "$nmatches topic(s)."
	} else {
		set curhelp(msg) "Not found."
	}
	seek $curhelp(fid) $startpos
	return $nmatches
}
