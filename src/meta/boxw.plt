# RCSid $Id: boxw.plt,v 1.1 2015/10/12 20:10:50 greg Exp $
#
#	Box and Whisker Plot.
#	This file represents a single curve with error
#	and extrema boxes.  The curve itself is A.
#	The upper box limit is B.  The lower box limit is C.
#	The maximum is D.  The minimum is E.
#

include = standard.plt		# standard definitions

fthick = 3			# frame thickness
othick = 0			# origin thickness
grid = 0			# grid on?

symfile = symbols.mta		# our symbol file

#
# Curve defaults:
#

Asymtype = "ex"		# symbol for main curve
Asymsize = 100 		# symbol radius
Athick = 2			# line thickness
Alintype = 1			# line type
Acolor = 4			# color

Bsymtype = "boxtop"		# upper box limit
Bsymsize = 150
Bthick = 0
Bcolor = 2

Csymtype = "boxbottom"		# lower box limit
Csymsize = 150
Cthick = 0
Ccolor = 2

Dsymtype = "uptee"		# maximum
Dsymsize = 150
Dthick = 0
Dcolor = 3

Esymtype = "downtee"		# minimum
Esymsize = 150
Ethick = 0
Ecolor = 3
