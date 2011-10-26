# RCSid: $Id: do_results.tcl,v 2.9 2011/10/26 20:06:17 greg Exp $
#
# Results screen for trad
#

set conv(TIFF-bw,nam)	"TIFF B&W"
set conv(TIFF-bw,com)	"ra_tiff -b %s %s"
set conv(TIFF-bw,suf)	.tif
set conv(TIFF-24,nam)	"TIFF 24-bit"
set conv(TIFF-24,com)	"ra_tiff %s %s"
set conv(TIFF-24,suf)	.tif
set conv(TIFF-Luv,nam)	"TIFF LogLuv"
set conv(TIFF-Luv,com)	"ra_tiff -L %s %s"
set conv(TIFF-Luv,suf)	.tif
set conv(GIF-bw,nam)	"GIF B&W"
set conv(GIF-bw,com)	"ra_gif -b %s %s"
set conv(GIF-bw,suf)	.gif
set conv(GIF-8,nam)	"GIF"
set conv(GIF-8,com)	"ra_gif -n 10 -d %s %s"
set conv(GIF-8,suf)	.gif
set conv(PPM-bin,nam)	"PPM (binary)"
set conv(PPM-bin,com)	"ra_ppm %s %s"
set conv(PPM-bin,suf)	.ppm
set conv(PPM-asc,nam)	"PPM (ASCII)"
set conv(PPM-asc,com)	"ra_ppm -a %s %s"
set conv(PPM-asc,suf)	.ppm
set conv(BMP-bw,nam)	"BMP B&W"
set conv(BMP-bw,com)	"ra_bmp -b %s %s"
set conv(BMP-bw,suf)	.bmp
set conv(BMP-clr,nam)	"BMP Color"
set conv(BMP-clr,com)	"ra_bmp %s %s"
set conv(BMP-clr,suf)	.bmp
set conv(PICT,nam)	"PICT 32-bit"
set conv(PICT,com)	"ra_pict -g 1.8 %s %s"
set conv(PICT,suf)	.pict
set conv(PS-bw,nam)	"PostScript B&W"
set conv(PS-bw,com)	"ra_ps -b -C %s %s"
set conv(PS-bw,suf)	.ps
set conv(PS-clr,nam)	"PostScript Color"
set conv(PS-clr,com)	"ra_ps -c -C %s %s"
set conv(PS-clr,suf)	.ps
set conv(tga-bw,nam)	"Targa B&W"
set conv(tga-bw,com)	"ra_t8 -b %s %s"
set conv(tga-bw,suf)	.tga
set conv(tga-8,nam)	"Targa 8-bit"
set conv(tga-8,com)	"ra_t8 -n 10 -d %s %s"
set conv(tga-8,suf)	.tga
set conv(tga-16,nam)	"Targa 16-bit"
set conv(tga-16,com)	"ra_t16 -2 %s %s"
set conv(tga-16,suf)	.tga
set conv(tga-24,nam)	"Targa 24-bit"
set conv(tga-24,com)	"ra_t16 -3 %s %s"
set conv(tga-24,suf)	.tga
set conv(types) {GIF-bw GIF-8 PICT PS-bw PS-clr BMP-bw BMP-clr PPM-asc\
		PPM-bin tga-bw tga-8 tga-16 tga-24 TIFF-bw TIFF-24}
set conv(typ) tga-24

proc testappend {flst tf} {	# test if tf exists and append to flst if so
	upvar $flst mylist
	if [file isfile $tf] {
		lappend mylist $tf
	}
}

proc list_views {} {		# List finished and unfinished pictures
	global radvar fvwbox ufvwbox alldone rawfroot
	set fpics {}
	set ufpics {}
	foreach vw $radvar(view) {
		if [file isfile ${rawfroot}_[lindex $vw 0].unf] {
			lappend ufpics [lindex $vw 0]
		} elseif {[file isfile $radvar(PICTURE)_[lindex $vw 0].hdr]} {
			lappend fpics [lindex $vw 0]
		}
	}
	$fvwbox delete 0 end
	eval $fvwbox insert end $fpics
	$ufvwbox delete 0 end
	eval $ufvwbox insert end $ufpics
	set alldone [expr [llength $fpics] == [llength $radvar(view)]]
}

proc delpic {} {		# Delete selected pictures
	global curmess
	set selected_pics [get_selpics 1]
	if {$selected_pics == {}} {
		set curmess "No pictures selected."
		return
	}
	if [tk_dialog .dlg {Verification} \
			"Really delete file(s) $selected_pics?" \
			questhead 0 {Delete} {Cancel}] {
		return
	}
	if {! [catch {eval exec rm $selected_pics < /dev/null} curmess]} {
		set curmess "Deleted [llength $selected_pics] file(s)."
	}
	list_views
}

proc get_selpics {{getall 0}} {		# return selected pictures
	global fvwbox ufvwbox radvar rawfroot
	set sl {}
	foreach i [$fvwbox curselection] {
		testappend sl $radvar(PICTURE)_[$fvwbox get $i].hdr
		if {$getall && $rawfroot != $radvar(PICTURE)} {
			testappend sl ${rawfroot}_[$fvwbox get $i].hdr
		}
		if {$getall && $radvar(ZFILE) != {}} {
			testappend sl $radvar(ZFILE)_[$fvwbox get $i].zbf
		}
	}
	foreach i [$ufvwbox curselection] {
		testappend sl ${rawfroot}_[$ufvwbox get $i].unf
		if {$getall && $radvar(ZFILE) != {}} {
			testappend sl $radvar(ZFILE)_[$ufvwbox get $i].zbf
		}
	}
	return $sl
}

proc dsppic {} {		# Display selected pictures
	global curmess dispcom radvar
	set selected_pics [get_selpics]
	if {$selected_pics == {}} {
		set curmess "No pictures selected."
		return
	}
	if {$radvar(EXPOSURE) == {}} {
		set ev 0
	} else {
		if [regexp {^[+-]} $radvar(EXPOSURE)] {
			set ev [expr {round($radvar(EXPOSURE))}]
		} else {
			set ev [expr {round(log($radvar(EXPOSURE))/log(2))}]
		}
	}
	foreach p $selected_pics {
		if {[string match *.unf $p] ||
				$radvar(PICTURE) == $radvar(RAWFILE)} {
			set dc [format $dispcom $ev $p]
		} else {
			set dc [format $dispcom 0 $p]
		}
		catch {eval exec $dc} curmess
	}
}

proc cnvpic {} {		# Convert selected pictures
	global curmess radvar conv convdest fvwbox
	set curmess "No finished pictures selected."
	foreach i [$fvwbox curselection] {
		set vw [$fvwbox get $i]
		set p $radvar(PICTURE)_$vw.hdr
		set df [format $convdest $vw]
		set curmess "Converting $p to $df..."
		update
		set cc [format $conv($conv(typ),com) $p $df]
		if {! [catch {eval exec $cc} curmess]} {
			set curmess "Done."
		}
	}
}

proc prtpic {} {		# Print selected pictures
	global curmess prntcom radvar fvwbox
	set curmess "No finished pictures selected."
	foreach i [$fvwbox curselection] {
		set p $radvar(PICTURE)_[$fvwbox get $i].hdr
		set curmess "Printing $p..."
		update
		set pc [format $prntcom $p]
		if {! [catch {eval exec $pc} curmess]} {
			set curmess "Done."
		}
	}
}

proc do_results w {		# Results screen
	global radvar curmess fvwbox ufvwbox dispcom prntcom conv \
			rawfroot convdest
	if {"$w" == "done"} {
		unset fvwbox ufvwbox convdest rawfroot
		return
	}
	frame $w
	# Finished view box
	label $w.vvl -text "Finished views"
	place $w.vvl -relx .0714 -rely .0610
	frame $w.vnl
	scrollbar $w.vnl.sb -relief sunken -command "$w.vnl.lb yview"
	listbox $w.vnl.lb -relief sunken -yscroll "$w.vnl.sb set" \
			-selectmode extended
	pack $w.vnl.sb -side right -fill y
	pack $w.vnl.lb -side left -expand yes -fill both
	place $w.vnl -relwidth .2143 -relheight .4268 -relx .0714 -rely .1220
	set fvwbox $w.vnl.lb
	helplink $w.vnl.lb trad results finished
	# Unfinished view box
	label $w.uvvl -text "Unfinished views"
	place $w.uvvl -relx .7143 -rely .0610
	frame $w.uvnl
	scrollbar $w.uvnl.sb -relief sunken -command "$w.uvnl.lb yview"
	listbox $w.uvnl.lb -relief sunken -yscroll "$w.uvnl.sb set" \
			-selectmode extended
	pack $w.uvnl.sb -side right -fill y
	pack $w.uvnl.lb -side left -expand yes -fill both
	place $w.uvnl -relwidth .2143 -relheight .4268 -relx .7143 -rely .1220
	set ufvwbox $w.uvnl.lb
	helplink $w.uvnl.lb trad results unfinished
	# Rescan button
	button $w.rsb -text Rescan -relief raised -command list_views
	place $w.rsb -relwidth .1071 -relheight .0610 -relx .4464 -rely .2439
	helplink $w.rsb trad results rescan
	# Delete button
	button $w.del -text Delete -relief raised -command delpic
	place $w.del -relwidth .1071 -relheight .0610 -relx .4464 -rely .3659
	helplink $w.del trad results delete
	# Display picture(s)
	if {! [info exists dispcom]} {
		set dispcom "ximage -e %+d %s >& /dev/null &"
	}
	button $w.dpyb -text Display -relief raised -command dsppic
	place $w.dpyb -relwidth .1071 -relheight .0610 -relx .0714 -rely .6098
	helplink $w.dpyb trad results display
	label $w.dpycl -text Command:
	place $w.dpycl -relx .2143 -rely .6098
	entry $w.dpyce -textvariable dispcom -relief sunken
	place $w.dpyce -relwidth .5714 -relheight .0610 -relx .3571 -rely .6098
	helplink $w.dpyce trad results dispcommand
	# Convert picture(s)
	button $w.cnvb -text Convert -relief raised -command cnvpic
	place $w.cnvb -relwidth .1071 -relheight .0610 -relx .0714 -rely .7317
	helplink $w.cnvb trad results convert
	menubutton $w.typb -text $conv($conv(typ),nam) -relief raised \
			-menu $w.typb.m
	place $w.typb -relwidth .1986 -relheight .0610 -relx .2143 -rely .7317
	helplink $w.typb trad results convtype
	menu $w.typb.m
	foreach t $conv(types) {
		$w.typb.m add radiobutton -variable conv(typ) -value $t \
			-label $conv($t,nam) \
			-command "$w.typb configure -text \"$conv($t,nam)\"
				set convdest $radvar(PICTURE)_%s$conv($t,suf)"
	}
	label $w.fil -text File:
	place $w.fil -relx .4486 -rely .7317
	set convdest $radvar(PICTURE)_%s$conv($conv(typ),suf)
	entry $w.file -textvariable convdest -relief sunken
	place $w.file -relwidth .4086 -relheight .0610 -relx .5200 -rely .7317
	helplink $w.file trad results convfile
	# Print picture(s)
	button $w.prtb -text Print -relief raised -command prtpic
	place $w.prtb -relwidth .1071 -relheight .0610 -relx .0714 -rely .8537
	helplink $w.prtb trad results print
	label $w.prtl -text Command:
	place $w.prtl -relx .2143 -rely .8537
	if {! [info exists prntcom]} {
		set prntcom "ra_ps %s | lpr"
	}
	entry $w.prte -textvariable prntcom -relief sunken
	place $w.prte -relwidth .5714 -relheight .0610 -relx .3571 -rely .8537
	helplink $w.prte trad results printcommand
	# Fill in views
	if {$radvar(RAWFILE) != {}} {
		set rawfroot $radvar(RAWFILE)
		if {$radvar(RAWFILE) == $radvar(PICTURE)} {
			set curmess "Warning: finished views are unfiltered"
		}
	} else {
		set rawfroot $radvar(PICTURE)
	}
	list_views
}
