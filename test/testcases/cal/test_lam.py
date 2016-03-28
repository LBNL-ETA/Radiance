# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class LamTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd):
		try:
			proc = self.call_one(cmd, 'call rlam', out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return lcompare.split_rad(raw)

	def test_lam(self):
		dat_de = ts.datafile('lam_de.dat')
		dat_en = ts.datafile('lam_en.dat')
		cmd = ['rlam', '-t:', dat_de, dat_en]
		result = self._runit(cmd)
		expect = [
		['eins:one'],
		['zwei:two'],
		['drei:three'],
		['vier:four'],
		['fuenf:five'],
	]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e),cmd))


# vi: set ts=4 sw=4 :
