
import os
import string

import unittest

def _find_raydir():
	try: thisfile = __file__
	except NameError:
		import support
		thisfile = support.__file__
	sdirn = string.count(thisfile, os.sep)
	if sdirn in [0, 1]:
		raydir = os.path.join('..','..')
	elif sdirn == 2:
		raydir = os.path.join('..')
	elif sdirn == 3:
		raydir = '.'
	else:
		unit_tools = os.path.split(thisfile)[0]
		py_tests = os.path.split(unit_tools)[0]
		test = os.path.split(py_tests)[0]
		raydir = os.path.split(test)[0]
	return raydir

RAYDIR = _find_raydir()
DATADIR = os.path.join(RAYDIR, 'test', 'test data')
BINDIR = os.path.join(RAYDIR, 'bin')

def binfile(fn):
	'''return the full/relative path to file in bin dir'''
	if os.name == 'nt':
		return os.path.normpath(os.path.join(BINDIR, fn + '.exe'))
	else: return os.path.normpath(os.path.join(BINDIR, fn))

def datafile(fn):
	'''return the full/relative path to file in data dir'''
	return os.path.normpath(os.path.join(DATADIR, fn))

def run_case(c):
	res = unittest.TestResult()
	s = unittest.makeSuite(c, 'test')
	s.run(res)
	if res.errors or res.failures:
		print ' failed'
		print '-----------------------------'
		for e in res.errors:
			print e[1]
		for e in res.failures:
			es = string.strip(e[1])
			sl = string.split(es, '\n')
			print sl[-1]
		print '-----------------------------'
	else:
		print ' ok'
		
