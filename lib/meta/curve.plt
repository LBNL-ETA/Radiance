#
#	This file defines a standard curve plot,
#	where a curve is represented by special
#	symbols at each point, with lines between.
#

include = standard.plt		# standard definitions

fthick = 3	# frame thickness
othick = 0	# origin thickness
lthick = 2	# line thickness
grid = 0	# grid on?

symfile = symbols.mta		# our symbol file
symsize = 100			# global symbol radius

#
# Curve defaults:
#

Asymtype = "triangle"		# symbol for A
Asymsize = symsize		# symbol radius for A
Athick = lthick			# line thickness for A
Alintype = 1			# line type for A
Acolor = 4			# color for A

Bsymtype = "square"
Bsymsize = symsize
Bthick = lthick
Blintype = 2
Bcolor = 2

Csymtype = "triangle2"
Csymsize = symsize
Cthick = lthick
Clintype = 3
Ccolor = 1

Dsymtype = "diamond"
Dsymsize = symsize
Dthick = lthick
Dlintype = 4
Dcolor = 3

Esymtype = "octagon"
Esymsize = symsize
Ethick = lthick
Elintype = 1
Ecolor = 2

Fsymtype = "crosssquare"
Fsymsize = symsize
Fthick = lthick
Flintype = 2
Fcolor = 4

Gsymtype = "exsquare"
Gsymsize = symsize
Gthick = lthick
Glintype = 3
Gcolor = 3

Hsymtype = "crossdiamond"
Hsymsize = symsize
Hthick = lthick
Hlintype = 4
Hcolor = 1
