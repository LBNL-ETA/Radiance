import os

from SCons.SConf import SConf # aka Configure

def find_radlib(env):
	v = env.FindFile('helvet.fnt', './lib')
	if not v:
		print '''
	Radiance auxiliary support files not found.
	-> Download from radiance-online.org and extract.
	'''
		env.Exit()

def find_x11(env):
	# Search for libX11, remember the X11 library and include dirs
	for d in ('/usr/X11R6', '/usr/X11', '/usr/openwin'):
		if os.path.isdir (d):
			incdir = os.path.join(d, 'include')
			libdir = os.path.join(d, 'lib')
			env.Prepend(CPPPATH=[incdir]) # add temporarily
			env.Prepend(LIBPATH=[libdir]) 
			conf = SConf(env)
			if conf.CheckLibWithHeader('X11', 'X11/X.h', 'C', autoadd=0):
				env.Replace(X11INCLUDE=incdir)
				env.Replace(X11LIB=libdir)
			env['CPPPATH'].remove(incdir) # not needed for now
			env['LIBPATH'].remove(libdir)
			if env['X11INCLUDE']:
				# Check for SGI stereo extension
				if conf.CheckCHeader('X11/extensions/SGIStereo.h'):
					env['RAD_STEREO'] = '-DSTEREO'
				else: env['RAD_STEREO'] = '-DNOSTEREO'
				env = conf.Finish ()
				break
			env = conf.Finish ()


def find_gl(env):
	# Check for libGL, set flag
	dl = [(None,None)] # standard search path
	if env.has_key('X11INCLUDE'): # sometimes found there (Darwin)
		dl.append((env['X11INCLUDE'], env['X11LIB']))
	for incdir, libdir in dl:
		if incdir: env.Prepend(CPPPATH=[incdir]) # add temporarily
		if libdir: env.Prepend(LIBPATH=[libdir])
		conf = SConf(env)
		if conf.CheckLibWithHeader('GL', 'GL/gl.h', 'C', autoadd=0):
			env['OGL'] = 1
		if incdir: env['CPPPATH'].remove(incdir) # not needed for now
		if libdir: env['LIBPATH'].remove(libdir)
		if env.has_key('OGL'):
			if incdir: env.Replace(OGLINCLUDE=[incdir])
			#if libdir: env.Replace(OGLLIB=[libdir])
			conf.Finish()
			break
		conf.Finish()


def find_libtiff(env):
	# Check for libtiff, set flag and include/lib directories
	dl = [ (None,None), ] # standard search path
	cfgi = env.get('TIFFINCLUDE')
	cfgl = env.get('TIFFLIB')
	if cfgi or cfgl:
		dl.insert(0,(cfgi, cfgl))
	for incdir, libdir in dl:
		if incdir: env.Prepend(CPPPATH=[incdir]) # add temporarily
		if libdir: env.Prepend(LIBPATH=[libdir])
		conf = SConf(env)
		if conf.CheckLib('tiff', 'TIFFInitSGILog',
				header='void TIFFInitSGILog(void);', autoadd=0):
			env['TIFFLIB_INSTALLED'] = 1
		if incdir: env['CPPPATH'].remove(incdir) # not needed for now
		if libdir: env['LIBPATH'].remove(libdir)
		if env.has_key('TIFFLIB_INSTALLED'):
			if incdir: env.Replace(RAD_TIFFINCLUDE=[incdir])
			if libdir: env.Replace(RAD_TIFFLIB=[libdir])
			conf.Finish()
			break
		conf.Finish()

