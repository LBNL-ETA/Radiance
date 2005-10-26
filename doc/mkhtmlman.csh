#!/bin/csh -fe
# RCSid $Id: mkhtmlman.csh,v 1.1 2005/10/26 23:56:21 greg Exp $
#
# Build man pages in HTML for Radiance website
#
cd man
mkdir man_html
foreach sect (1 3 5)
cd man$sect
foreach mpage (*.$sect)
	nroff -man $mpage \
		| man2html -title "Radiance $mpage:r program" -nodepage \
			     -cgiurl '$title.$section.html' \
		> ../man_html/$mpage.html
end
cd ..
end
tar czf ../man_html.tar.gz man_html
echo "Pages in man_html.tar.gz"
exec rm -r man_html
