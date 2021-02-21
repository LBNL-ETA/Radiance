#
#	This file defines a standard scatter plot,
#	where a curve is represented by special
#	symbols at each point.
#

include = standard.plt		# standard definitions

symfile = symbols.mta		# our symbol file
symsize = 100			# global symbol radius

#
# Curve defaults:
#

Asymtype = "triangle"		# symbol for A
Asymsize = symsize		# symbol radius for A
Athick = 0			# line thickness for A
Acolor = 4			# color for A

Bsymtype = "square"
Bsymsize = symsize
Bthick = 0
Bcolor = 2

Csymtype = "triangle2"
Csymsize = symsize
Cthick = 0
Ccolor = 1

Dsymtype = "diamond"
Dsymsize = symsize
Dthick = 0
Dcolor = 3

Esymtype = "octagon"
Esymsize = symsize
Ethick = 0
Ecolor = 4

Fsymtype = "crosssquare"
Fsymsize = symsize
Fthick = 0
Fcolor = 2

Gsymtype = "exsquare"
Gsymsize = symsize
Gthick = 0
Gcolor = 1

Hsymtype = "crossdiamond"
Hsymsize = symsize
Hthick = 0
Hcolor = 3
