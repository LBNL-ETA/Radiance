#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' run_tests.py - Run Radiance test suite
2013 - 2016 Georg Mischler

Invocation:
  As script, call with -H for instructions.
  As module, see docstrings in class RadianceTests for instructions.
'''
from __future__ import division, print_function, unicode_literals
__all__ = ['TESTCATS', 'RadianceTests', 'main']
import os
import sys
import types
import argparse
import unittest

SHORTPROGN = os.path.splitext(os.path.basename(sys.argv[0]))[0]
TESTCATS = ('pyrad','cal','cv','gen','hd','meta','ot','px', 'rt','util',)

class Error(Exception): pass

class RadianceTests():
	'''Test Radiance programs below subdirectory "testcases" recursively.

	This class will create a virtual module "testsupport" with constants
	and functions for use in individual test cases.

	  bindir   - list of paths added to PATH
	  path     - complete list of paths currently in PATH
	  radlib   - list of paths added to RAYPATH
	  raypath  - complete list of paths currently in RAYPATH
	  pyradlib - absolute path of python support library directory
	  datadir  - absolute path of test data directory
	  datafile([sub, [sub,]] fn) - return absolute file path in
		data directory or in a subdirectory
	'''
	def __init__(self, **args):
		'''Args are:
		   bindir=[directory ...] - will be prepended to PATH for tests
		   radlib=[directory ...] - will be prepended to RAYPATH for tests
		   cat=[category ...]     - only test those categories (else TESTCATS)
		   V=False                - if True, verbose listing of executed tests
		'''
		self.bindir = args.get('bindir')
		self.radlib = args.get('radlib')
		self.testcats = args.get('cat')
		self.V = 2 if args.get('V') else 1
		try:
			old_osenv = os.environ

			thisfile = os.path.abspath(__file__)
			self.thisdir = os.path.split(thisfile)[0]
			self.testdir = os.path.join(self.thisdir, 'testcases')
			old_syspath = sys.path
			if self.testdir not in sys.path:
				sys.path.insert(0, self.testdir)

			oldpath = old_osenv.get('PATH', '')
			pathlist = oldpath.split(os.pathsep)
			if self.bindir:
				for bdir in self.bindir:
					if bdir not in pathlist:
						pathlist.insert(0, bdir)
				os.environ['PATH'] = os.pathsep.join(pathlist)

			oldraypath = old_osenv.get('RAYPATH', '')
			raypathlist = oldraypath.split(os.pathsep)
			if self.radlib:
				for rlib in self.radlib:
					if rlib not in raypathlist:
						raypathlist.insert(0, rlib)
				os.environ['RAYPATH'] = os.pathsep.join(raypathlist)

			pyradlib = None
			for rp in raypathlist:
				prd = os.path.join(rp, 'pyradlib')
				#print(prd)
				if os.path.isdir(prd) or os.path.islink(prd):
					#print('--- found !!!')
					if rp not in sys.path:
						sys.path.insert(0, rp)
					pyradlib = prd
			if not pyradlib:
				raise Error('Python support library directory "pyradlib" not found on RAYPATH')

			# Since this here is the best place to figure out where
			# everything is, we create an ad-hoc module to hold directory
			# paths and functions for creating file paths.
			# The individual test cases can then "import testsupport" to
			# get access to the complete information.
			supmod = types.ModuleType(str('testsupport'))
			datadir = os.path.join(self.thisdir, 'test data')
			supmod.bindir = self.bindir
			supmod.datadir = datadir
			supmod.datafile = lambda fn,*ffn: os.path.join(datadir, fn, *ffn)
			supmod.path = pathlist
			supmod.radlib = self.radlib
			supmod.raypath = raypathlist
			supmod.pyradlib = pyradlib
			sys.modules['testsupport'] = supmod

			self._run_all()
		finally:
			if hasattr(sys.modules, 'testsupport'):
				del sys.modules['testsupport']
			os.environ = old_osenv
			sys.path = old_syspath

	def _run_all(self):
		runner = unittest.TextTestRunner(verbosity=self.V)
		if self.testcats:
			cats = self.testcats
		else: cats = TESTCATS
		for dir in cats:
			loader = unittest.defaultTestLoader
			fulldir = os.path.join(self.testdir, dir)
			if not os.path.isdir(fulldir):
				raise Error('No such directory: "%s"' % fulldir)
			suite = loader.discover(fulldir, 'test*.py',
					top_level_dir=self.thisdir)
			count = suite.countTestCases()
			if count:
				print('\n--', dir, '----', file=sys.stderr)
				runner.run(suite)


def main():
	'''Main function for invocation as script. See usage instructions with -H'''
	parser = argparse.ArgumentParser(add_help=False,
		description='Run Radiance test suite',)
	parser.add_argument('-V', action='store_true',
		help='Verbose: Print all executed test cases to stderr')
	parser.add_argument('-H', action='help',
		help='Help: print this text to stderr and exit')
	parser.add_argument('-p', action='store', nargs=1,
		dest='bindir', metavar='bindir', help='Path to Radiance binaries')
	parser.add_argument('-l', action='store', nargs=1,
		dest='radlib', metavar='radlib', help='Path to Radiance library')
	parser.add_argument('-c', action='store', nargs=1,
		dest='cat', metavar='cat', help='Category of tests to run (else all)')
	args = parser.parse_args()
	RadianceTests(**vars(args))


if __name__ == '__main__':
	try: main()
	except KeyboardInterrupt:
		sys.stderr.write('*cancelled*\n')
		exit(1)
	except Error as e:
		sys.stderr.write('%s: %s\n' % (SHORTPROGN, str(e)))
		exit(-1)

# vi: set ts=4 sw=4 :
