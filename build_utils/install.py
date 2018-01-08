from __future__ import division, print_function, unicode_literals

import os
import stat

def install_dir(env, instlist, sbase, tbase, subdir):
	fulldir = os.path.join(sbase, subdir)
	l = os.listdir(fulldir)
	instdir = os.path.join(tbase, subdir)
	for item in l:
		if (item[0] == '.'
			or item[-1] in '%~'
			or item == 'CVS'
			or item.lower().startswith('cmake')):
			continue
		elif os.path.isdir(os.path.join(fulldir, item)):
			install_dir(env, instlist, sbase, tbase, os.path.join(subdir, item))
		else:
			inst = env.Install(instdir, os.path.join(fulldir, item))
			instlist.append(inst)

def install_rlibfiles(env):
	buildrlib = env['RAD_BUILDRLIB']
	if buildrlib[0] == '#':
		buildrlib = buildrlib[1:]
	sbase = os.path.join(os.getcwd(), buildrlib)
	instlist = []
	install_dir(env, instlist, sbase, env['RAD_RLIBDIR'], '')
	env.Append(RAD_RLIBINSTALL=instlist)


def install_manfiles(env):
	buildman = env['RAD_BUILDMAN']
	if buildman[0] == '#':
		buildman = buildman[1:]
	sbase = os.path.join(os.getcwd(), buildman)
	il = []
	install_dir(env, il, sbase, env['RAD_MANDIR'], '')
	env.Append(RAD_MANINSTALL=il)


def install_script(target, source, env):
	#for t,s in map(None,target,source):
	for t,s in zip(target, source):
		sf = open(str(s), 'r')
		tf = open(str(t), 'w')
		tf.write(sf.read())
		os.chmod(str(t),
				stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
				stat.S_IRGRP | stat.S_IWGRP |
				stat.S_IROTH | stat.S_IWOTH)

def install_tclscript(target, source, env):
	# XXX posix only for the moment
	#for t,s in map(None, target,source):
	for t,s in zip(target, source):
		bindir = os.path.split(t.get_abspath())[0]
		instdir = os.path.split(bindir)[0]
		tcllibdir = os.path.join(instdir, 'lib', 'tcl')
		sf = open(str(s), 'r')
		tf = open(str(t), 'w')
		# as recommended in the wish manpage
		tf.write('''#!/bin/sh\n''')
		for line in sf.readlines():
			if line.startswith('#!'):
				tf.write('### ' + line)
			if line.startswith('set radlib'):
				# XXX this should really be handled by an envvar.
				line = 'set radlib %s\n' % tcllibdir
			tf.write(line)
		os.chmod(str(t), 
				stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
				stat.S_IRGRP | stat.S_IWGRP |
				stat.S_IROTH | stat.S_IWOTH)

def build_with_pyinstaller(targets, sources, env):
	workpath = env.Dir('$RAD_BUILDOBJ', 'pybuild') # --workpath @
	specpath = workpath # --specpath @
	# -F   # one-file
	# --noconsole
	pass

# vi: set ts=4 sw=4 :
