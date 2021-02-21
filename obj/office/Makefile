#
#  A basic office scene
#

VIEW = -vp 8 36 -27 -vd -.46 -.23 .79 -vu 0 1 0 -vh 39.9 -vv 27.5
DEV = x11

view1:	model.oct
	rview $(VIEW) -o $(DEV) -av .1 .1 .1 model.oct

display:	/usr2/greg/pic/model.new
	ximage =+50+50 -f /usr2/greg/pic/model.new

model.oct:	modelb.oct window blinds lights lamp picture
	oconv -f -i modelb.oct window blinds lights lamp picture > model.oct

modelb.oct:	model.b90 desk misc bookbox.norm chair1.norm chair2.norm \
cornerdesk.norm mug.norm openbook.norm pageverts
	oconv model.b90 desk misc > modelb.oct
