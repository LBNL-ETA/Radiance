# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class GetinfoTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd):
		try:
			proc = self.call_one(cmd, 'call getinfo', out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return lcompare.split_headers(raw)

	def test_getinfo(self):
		picfile = ts.datafile('Earth128.pic')
		cmd = ['getinfo', picfile]
		result = self._runit(cmd)
		exps = '''%s:
	Xim format conversion by:
	FORMAT=32-bit_rle_rgbe
	pfilt -e 2 -x 512 -y 512 -p 1 -r .67
	EXPOSURE=4.926198e+00
	normpat
	pfilt -1 -e .2
	EXPOSURE=2.000000e-01
	pfilt -x 128 -y 128
	PIXASPECT=0.500000
	EXPOSURE=2.571646e+00
	
''' % picfile
		expect = lcompare.split_headers(exps)
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_getinfo_d(self):
		picfile = ts.datafile('Earth128.pic')
		cmd = ['getinfo', '-d', picfile]
		result = self._runit(cmd)
		exps = '''%s: -Y 128 +X 128\n''' % picfile
		expect = lcompare.split_headers(exps)
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e),cmd))


# vi: set ts=4 sw=4 :
