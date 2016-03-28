# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class XformTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd):
		try:
			proc = self.call_one(cmd, 'call xform', out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return lcompare.split_rad(raw)

	def test_xform_e_s(self):
		cmd = 'xform -e -s 3.14159'.split() + [ts.datafile('xform_1.dat')]
		result = self._runit(cmd)
		expect = lcompare.split_radfile(ts.datafile('xform_res1.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_xform_e_mx(self):
		cmd = 'xform -e -mx'.split() + [ts.datafile('xform_2.dat')]
		result = self._runit(cmd)
		expect = lcompare.split_radfile(ts.datafile('xform_res2.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e),cmd))


# vi: set ts=4 sw=4 :
