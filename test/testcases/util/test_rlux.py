# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class RluxTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd, data):
		failmsg = None
		try:
			proc = self.call_one(cmd, 'test rlux',
					_in=PIPE, out=PIPE, universal_newlines=True)
			raw, err = proc.communicate(data)
		except Error as e:
			failmsg = str(e)
		if failmsg:
			self.fail(failmsg)
		return lcompare.split_rad(raw)

	def test_rlux(self):
		octfn = ts.datafile('corner', 'corner.oct')
		data = '''0.5 0.5 1.5 0 0 -1
			0.3 0.7 0.5 0 0 -1
			0.7 0.3 0.5 0 0 -1
			0.5 0.5 0.5 0 1 0
			0.5 0.5 0.5 -1 0 0
			-1 0.5 0.5 1 0 0
			0.5 2 0.5 0 -1 0
				'''
		exp = [[37796.5481], [0], [37796.5481], [0], [0], [46062.1953], [0], []]
		cmd = ['rlux', octfn]
		resl = self._runit(cmd, data)
		lcompare.llcompare(resl, exp)



# vi: set ts=4 sw=4 :
