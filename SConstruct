
import os
import sys
import string

OPTFILE = 'rayopts.py'
SourceSignatures('timestamp')
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

if os.name == 'posix':
	from build_utils import install
	csh_b = Builder(action = install.install_cshscript,
			suffix = '', src_suffix = '.csh')
	env.Append(BUILDERS={'InstallCsh': csh_b})

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

Default('.')

# virtual targets
env.Alias('bininstall',  '$RAD_BINDIR')
env.Alias('rlibinstall', '$RAD_RLIBDIR')
env.Alias('maninstall',  '$RAD_MANDIR')

env.Alias('build',   ['#bin'])
env.Alias('test',    ['#src/test'])
env.Alias('install', ['bininstall', 'rlibinstall', 'maninstall'])

