
import os
import sys
import string
import ConfigParser


_platdir = 'platform'


def POSIX_setup(env):
    # common stuff for all posix systems
    env['RAD_PROCESS'] = string.split('unix_process.c')


def read_plat(env, args, fn):
	cfig = ConfigParser.ConfigParser(env.Dictionary())
	cfig.read(fn)
	buildvars = [['CC'], # replace
			['CPPPATH', 'CPPDEFINES', 'CPPFLAGS', 'CCFLAGS',
			'LIBPATH', 'LINKFLAGS']] # append
	vars = [
		['install',
			['RAD_BASEDIR', 'RAD_BINDIR', 'RAD_RLIBDIR', 'RAD_MANDIR'],
			[]],
		['code',
			['RAD_SPEED'],
			['RAD_COMPAT', 'RAD_MEMCOMPAT', 'RAD_MATHCOMPAT', 'RAD_ARGSCOMPAT',
			'RAD_MLIB', 'RAD_PROCESS']],
	]
	if args.get('RAD_DEBUG',0):
		vars.insert(0, ['debug'] + buildvars)
	else: vars.insert(0, ['build'] + buildvars)
	for section in vars:
		if cfig.has_section(section[0]):
			for p in section[1]:
				try: v = cfig.get(section[0], p)
				except ConfigParser.NoOptionError: continue
				env[p] = v
				#print '%s: %s' % (p, env[p])
			for p in section[2]:
				try: v = cfig.get(section[0], p)
				except ConfigParser.NoOptionError: continue
				apply(env.Append,[],{p:string.split(v)})
				#print '%s: %s' % (p, env[p])
	# XXX Check that basedir exists.
	for k in ['RAD_BINDIR', 'RAD_RLIBDIR', 'RAD_MANDIR']:
		if (env.has_key('RAD_BASEDIR') and env.has_key(k)
				and not os.path.isabs(env[k])):
			env[k] = os.path.join(env['RAD_BASEDIR'],env[k])


def load_plat(env, args, platform=None):
	if os.name == 'posix':
		POSIX_setup(env)
	if platform == None: # override
		p = sys.platform
	else: p = platform
	pl = []
	print 'Detected platform "%s" (%s).' % (sys.platform, os.name)
	for i in [len(p), -1, -2]:
		pfn = os.path.join(_platdir, p[:i] + '_custom.cfg')
		if os.path.isfile(pfn):
			print 'Reading configuration "%s"' % pfn
			read_plat(env, args, pfn)
			return 1
		pfn = os.path.join(_platdir, p[:i] + '.cfg')
		if os.path.isfile(pfn):
			print 'Reading configuration "%s"' % pfn
			read_plat(env, args, pfn)
			return 1
		
	if os.name == 'posix':
		pfn = os.path.join(_platdir, 'posix.cfg')
		if os.path.isfile(pfn):
			print 'Reading generic configuration "%s".' % pfn
			read_plat(env, args, pfn)
			return 1

	print 'Platform "%s/%s" not supported yet' % (os.name, sys.platform)
	sys.exit(2)

	
	
