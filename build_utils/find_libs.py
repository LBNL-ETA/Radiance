import os

from SCons.SConf import SConf # aka Configure

def find_x11(env):
	# Search for libX11, remember the X11 library and include dirs
	for d in ('/usr/X11R6', '/usr/X11', '/usr/openwin'):
		if os.path.isdir (d):
			incdir = os.path.join(d, 'include')
			libdir = os.path.join(d, 'lib')
			env.Append(CPPPATH=[incdir]) # add temporarily
			env.Append(LIBPATH=[libdir]) 
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
	conf = SConf(env)
	if conf.CheckLibWithHeader('GL', 'GL/gl.h', 'C', autoadd=0):
		env['OGL'] = 1
	elif env.has_key('X11INCLUDE'): # sometimes found there
		incdir = env['X11INCLUDE']
		libdir = env['X11LIB']
		env.Append(CPPPATH=[incdir]) # add temporarily
		#env.Append(LIBPATH=[libdir])
		if conf.CheckLibWithHeader('GL', 'GL/gl.h', 'C', autoadd=0):
			env['OGL'] = 1
			env.Replace(OGLINCLUDE=incdir)
			env.Replace(OGLLIB=libdir)
		if incdir: env['CPPPATH'].remove(incdir) # not needed for now
		#if libdir: env['LIBPATH'].remove(libdir)
	conf.Finish()



