# RCSid: $Id: util.tcl,v 2.4 2003/02/22 02:07:30 greg Exp $
#
# General utility routines
#

proc beep {} {		# ring the bell
	puts -nonewline stderr "\a"
}

proc view_txt fn {		# view a text file
	global env
	if [info exists env(EDITOR)] {
		if {[lsearch -exact {vi vedit view ex edit ed red} \
				[file tail [lindex $env(EDITOR) 0]]] >= 0} {
			eval exec xterm +sb -e $env(EDITOR) $fn
		} else {
			eval exec $env(EDITOR) $fn
		}
	} else {
		foreach f $fn {
			exec xedit $f
		}
	}
}

proc cardval {var cvl} {	# set variable to one of a set of values
	upvar $var v
	set tv [string toupper $v]
	foreach cv $cvl {
		if {[string first $tv [string toupper $cv]] == 0} {
			set v $cv
			return $cv
		}
	}
	return {} ;		# this would seem to be an error
}

proc writevars {f vl} {		# write variables to a file
	foreach v $vl {
		upvar $v myv
		if [catch {set il [array names myv]}] {
			puts $f "set $v [list $myv]"
		} else {
			foreach i $il {
				puts $f "set ${v}($i) [list $myv($i)]"
			}
		}
	}
}
