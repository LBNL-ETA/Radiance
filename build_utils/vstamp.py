
import getpass
import socket
import string
import time

# create a version stamp
def build_version_c(target, source, env):
	s = open (str(source[0]), 'r')
	radver = string.strip (s.readline ())
	s.close()
	date = time.ctime (time.time())
	user = getpass.getuser()
	sysname = socket.gethostname()
	t = open (str(target[0]), 'w')
	t.write (
		'char VersionID[]="%s lastmod %s by %s on %s";\n' %
		(radver, date, user, sysname))
	t.close()

