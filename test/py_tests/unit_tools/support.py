
import os
import sys
import string

import unittest

def _find_raydir():
	if __name__ == '__main__':
		thisfile = sys.argv[0]
	else: thisfile = sys.modules[__name__].__file__
	unit_tools = os.path.abspath(os.path.dirname(thisfile))
	py_tests = os.path.dirname(unit_tools)
	test = os.path.dirname(py_tests)
	raydir = os.path.dirname(test)
	return raydir

RAYDIR = _find_raydir()
DATADIR = os.path.join(RAYDIR, 'test', 'test data')
BINDIR = os.path.join(RAYDIR, 'bin')

def binfile(fn):
	'''return the full path to file in bin dir'''
	if os.name == 'nt':
		return os.path.normpath(os.path.join(BINDIR, fn + '.exe'))
	else: return os.path.normpath(os.path.join(BINDIR, fn))

def datafile(fn):
	'''return the full path to file in data dir'''
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
		
