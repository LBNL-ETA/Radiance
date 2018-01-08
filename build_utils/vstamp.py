from __future__ import division, print_function, unicode_literals

import getpass
import socket
import time

# create a version stamp
def build_version_c(target, source, env):
	s = open (str(source[0]), 'r')
	radver = s.readline().strip()
	s.close()
	date = time.ctime (time.time())
	user = getpass.getuser()
	sysname = socket.gethostname()
	t = open (str(target[0]), 'w')
	t.write (
		'char VersionID[]="%s lastmod %s by %s on %s";\n' %
		(radver, date, user, sysname))
	t.close()

# vi: set ts=4 sw=4 :

