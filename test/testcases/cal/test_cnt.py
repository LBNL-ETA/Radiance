# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import itertools
import unittest

from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class CntTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd, actstr):
		try:
			proc = self.call_one(cmd, actstr, out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return [s.strip() for s in raw.split('\n') if s]

	def test_cnt_1d(self):
		for n in (5,10,100,1000,5000, 10000, 99999):
			cmd = ['cnt', str(n)]
			exp = range(n)
			raw = self._runit(cmd, 'test cnt 1dim')
			res = [float(s) for s in raw if s]
			try: lcompare.lcompare(res, exp)
			except lcompare.error as e:
				self.fail(('test_cnt_1 n=%d - ' % n) + str(e))

	def test_cnt_2d(self):
		# values higher than 1000/1000 will take rather long
		for n,m in ((3,3), (10,5), (100,200), ):#(1000,200)):
			cmd = ['cnt', str(n), str(m)]
			exp = itertools.product(range(n), range(m))
			raw = self._runit(cmd, 'test cnt 2dim')
			res = [[float(ss) for ss in s.split('\t')] for s in raw]
			try: lcompare.llcompare(res, exp)
			except lcompare.error as e:
				self.fail(('test_cnt_2 %d/%d - ' % (n,m)) + str(e))

	def test_cnt_3d(self):
		# values higher than 100/100/100 will take rather long
		for n,m,o in ((3,3,3), (10,5,7), (44, 33,22),):# (200,50,10)):
			cmd = ['cnt', str(n), str(m), str(o)]
			exp = itertools.product(range(n), range(m), range(o))
			res = self._runit(cmd, 'test cnt 3dim')
			res = [[float(ss) for ss in s.split('\t')] for s in res]
			try: lcompare.llcompare(res, exp)
			except lcompare.error as e:
				self.fail(('test_cnt_3 %d/%d/%d - ' % (n,m,o)) + str(e))


# vi: set ts=4 sw=4 :
