
import os
import sys
import string

OPTFILE = 'rayopts.py'

def set_opts(env):
    # XXX add some caching
    opts = Options(OPTFILE, ARGUMENTS)
    opts.Add('SKIP', 'Skip Display of License terms', 0)
    opts.Add('RAD_BINDIR',  'Install executables here',   env['RAD_BINDIR'])
    opts.Add('RAD_RLIBDIR', 'Install support files here', env['RAD_RLIBDIR'])
    opts.Add('RAD_MANDIR',  'Install man pages here',     env['RAD_MANDIR'])
    opts.Add('RAD_DEBUG',   'Build a debug version',  0)
    opts.Update(env) 
    opts.Save(OPTFILE, env)
    Help(opts.GenerateHelpText(env, sort=cmp))
	# where stuff is located in the source tree
    env['RAD_BUILDLIB']  = '#src/lib'
    env['RAD_BUILDBIN']  = '#bin'
    env['RAD_BUILDRLIB'] = '#lib'
    env['RAD_BUILDMAN']  = '#doc/man'

def allplats_setup(env):
    from build_utils import find_libs
    find_libs.find_x11(env)
    find_libs.find_gl(env) # OpenGL
    #find_libs.find_pixar(env) # PIXAR_LIB for src/px/ra_pixar.c

def post_common_setup(env):
    env.Append(CPPPATH = [os.path.join('#src', 'common')])
    env.Append(LIBPATH=['../lib']) # our own libs
    if not env.has_key('RAD_MLIB'):
        env['RAD_MLIB'] = [] #['m'] # no seperate mlib on Win

def shareinstall_setup(env):
	from build_utils import install
	# only scan for those files when we actually need them
	if 'install' in sys.argv or 'rlibinstall' in sys.argv:
		install.install_rlibfiles(env)
	if 'install' in sys.argv or 'maninstall' in sys.argv:
		install.install_manfiles(env)

# Set up build environment
env = Environment()

# configure platform-specific stuff
from build_utils import load_plat
load_plat.load_plat(env, ARGUMENTS, platform=None)

# override options
set_opts(env)

# accept license
if not env['SKIP'] and not '-c' in sys.argv:
    from build_utils import copyright
    copyright.show_license()

# fill in generic config
allplats_setup(env)

# Bring in all the actual things to build
Export('env')
SConscript(os.path.join('src', 'common', 'SConscript'))
post_common_setup(env)
for d in Split('meta cv gen ot rt px hd util cal'):
    print d
    SConscript(os.path.join('src', d, 'SConscript'))

if string.find(string.join(sys.argv[1:]), 'install') > -1:
	shareinstall_setup(env)

# virtual targets
# RAD_XXXINSTALL are filled in by the local scripts
env.Alias('bininstall',  env.get('RAD_BININSTALL', []))
env.Alias('rlibinstall', env.get('RAD_RLIBINSTALL',[]))
env.Alias('maninstall',  env.get('RAD_MANINSTALL', []))

env.Alias('build',   ['#bin'])
env.Alias('test',    ['#src/test'])
env.Alias('install', ['bininstall', 'rlibinstall', 'maninstall'])

# Further virtual targets are defined locally:
# meta_special: mt1601 okimate imagew mt160 mx80 impress aed5
#               tcurve tscat tbar mtext libt4014.a plotout t4014
# px_special:   ra_im, t4027, paintjet, mt160t, greyscale, colorscale, d48c
# util_special: scanner, makedist (not for Windows yet)
env.Alias('special', ['meta_special', 'px_special', 'util_special'])
env.Alias('special_install', ['meta_special_install',
		'px_special_install', 'util_special_install'])

