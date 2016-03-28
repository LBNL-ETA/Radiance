# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import math
import unittest

from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class EvTestCase(unittest.TestCase, ProcMixin):

	ltest = [
		['2.3 + 5 * 21.7 / 1.43', [2.3 + 5 * 21.7 / 1.43]],
		['if(1, 1, 0)', [1]],
		['if(0, 1, 0)', [0]],
		['select(3, 1, 2, 3, 4, 5)', [3]],
		['floor(5.743)', [5]],
		['ceil(5.743)', [6]],
		['sqrt(7.4)', [math.sqrt(7.4)]],
		['exp(3.4)', [math.exp(3.4)]],
		['log(2.4)', [math.log(2.4)]],
		['log10(5.4)', [math.log10(5.4)]],
		['sin(.51)', [math.sin(.51)]],
		['cos(.41)', [math.cos(.41)]],
		['tan(0.77)', [math.tan(0.77)]],
		['asin(0.83)', [math.asin(0.83)]],
		['acos(.94)', [math.acos(.94)]],
		['atan(0.22)', [math.atan(0.22)]],
		['atan2(0.72, 0.54)', [math.atan2(0.72, 0.54)]],
	]

	def _runit(self, cmd):
		try:
			proc = self.call_one(cmd, 'call ev', out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return raw.strip()

	def test_ev_singleres(self):
		for expr, expect in self.ltest:
			cmd = ['ev', expr]
			result = self._runit(cmd)
			try: lcompare.lcompare([result], expect)
			except lcompare.error as e:
				self.fail('%s [%s]' % (str(e),cmd))

	def x_test_ev_multipleres(self):
		pass # XXX implement


# vi: set ts=4 sw=4 :
