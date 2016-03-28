# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

	# We currently only test for output of recno after receiving a
	# single column list from cnt().
	# This serves to uncover the nl eating bug in buffered pipes
	# in Visual Studio 2015 before update 2.

import struct
import unittest

from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class RcalcTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd1, cmd2, actstr, nl=False):
		try:
			procs = self.call_many((cmd1, cmd2), actstr, out=PIPE,
					universal_newlines=nl)
			res, eres = procs[-1].communicate()
			if eres: print('Errors', eres) # XXX
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			for proc in procs:
				proc.wait()
		if nl: return [int(s) for s in res.split()]
		return res
		
		#return lcompare.split_rad(raw)

	def test_recno_oascii(self):
		for n in (5, 10, 55, 200,1321, 1328,1329,1330,1331,1332, 2000, 3721, 5000,9876):
			orig = range(n)
			exp = [i+1 for i in orig]
			cntcmd = ['cnt', str(n)]
			aacmd = ['rcalc', '-e', '$1=recno']
			aresl = self._runit(cntcmd, aacmd, 'calculate recno', nl=True)
			try: lcompare.lcompare(aresl, exp)
			except lcompare.error as e:
				self.fail(('recno_oascii n=%d -- ' % n) +str(e))

	def test_recno_ofloat(self):
		for n in (5, 10, 55, 200,1321, 1328,1329,1330,1331,1332, 2000, 3721, 5000):
			orig = range(n)
			exp = (i+1 for i in orig)
			cntcmd = ['cnt', str(n)]
			afcmd = ['rcalc', '-of', '-e', '$1=recno']
			res = self._runit(cntcmd, afcmd, 'calculate recno')
			if len(res) != n * 4:
				self.fail(('recno_ofloat n=%d -- Length of resulting data '
							'differs from expected  (%d != %d)')
						% (n, len(res), n*4))
			iresl = struct.unpack('f'*n, res)
			try: lcompare.lcompare(iresl, exp)
			except lcompare.error as e:
				self.fail(('recno_ofloat n=%d -- ' % n) +str(e))

	def test_recno_odouble(self):
		for n in (5, 10, 55, 200,1321, 1328,1329,1330,1331,1332, 2000, 3721, 5000):
			orig = range(n)
			exp = (i+1 for i in orig)
			cntcmd = ['cnt', str(n)]
			adcmd = ['rcalc', '-od', '-e', '$1=recno']
			res = self._runit(cntcmd, adcmd, 'calculate recno')
			if len(res) != n * 8:
				self.fail(('recno_odouble n=%d -- Length of resulting data '
							'differs from expected  (%d != %d)')
						% (n, len(res), n*8))
			iresl = struct.unpack('d'*n, res)
			try: lcompare.lcompare(iresl, exp)
			except lcompare.error as e:
				self.fail(('recno_odouble n=%d -- ' % n) +str(e))


# vi: set ts=4 sw=4 :
