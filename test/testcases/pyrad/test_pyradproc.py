# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import struct
import operator
import unittest
from functools import reduce

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class PyradprocTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmdl, actstrl, _in=None, indata=None, nl=False):
		try:
			procs = []
			if isinstance(cmdl[0], (type(b''), type(u''))):
				proc = self.call_one(cmdl, actstrl, _in=_in, out=PIPE,
						universal_newlines=nl)
				procs = [proc]
			elif len(cmdl) == 2:
				procs = self.call_two(cmdl[0], cmdl[1], actstrl[0], actstrl[1],
						_in=_in, out=PIPE, universal_newlines=nl)
			else:
				procs = self.call_many(cmdl, actstrl, _in=_in, out=PIPE,
						universal_newlines=nl)
			if indata: # only for small amounts of data
				procs[0].stdin.write(indata)
				procs[0].stdin.close()
				res = procs[-1].stdout.read()
			else:
				res, eres = procs[-1].communicate()
				if eres: print('Errors', eres) # XXX
		except Error as e:
			raise
			self.fail('%s [%s]' % (str(e), self.qjoin(cmdl)))
		finally:
			for proc in procs:
				proc.wait()
		if nl: return filter(None, lcompare.split_rad(res))
		return res
		

	def test_call_1_text(self):
		for n in (5, 10, 55, 200,1328,1329,1330,1331,1332, 2000, 5000,99999):
			exp = range(n)
			cntcmd = ['cnt', str(n)]
			res = self._runit(cntcmd, 'run one process', nl=True)
			resl = [r[0] for r in res]
			try: lcompare.lcompare(resl, exp)
			except lcompare.error as e:
				print(resl, exp)
				self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_call_1_bin(self):
		exp = (2160, 8640, 4320, 1080, 7560)
		cntcmd = ['total', '-of']
		tfn = ts.datafile('histo.dat')
		res = self._runit(cntcmd, 'run one process', _in=tfn, nl=False)
		resl = struct.unpack('f'*len(exp), res)
		try: lcompare.lcompare(resl, exp)
		except lcompare.error as e:
				self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_call_1_noexe(self):
		# let's hope a program with that name doesn't exist
		cmd = ['iswaxubdzkcgrhmltqvyfepojn', '-x', '-y', '-z']
		self.assertRaises(Error, self._runit, cmd, 'run one process')

	def test_call_2_text(self):
		# with values higher than 44721, total will return a float in e-format.
		for n in (5, 10, 55, 200,1328,1329,1330,1331,1332, 2000, 5000,
				44721, 44722, 99999):
			exp = [reduce(operator.add, range(n))]
			cntcmd = ['cnt', str(n)]
			totcmd = ['total']
			resl = self._runit([cntcmd, totcmd], ['run p1','run p2'], nl=True)
			try: lcompare.lcompare(resl, [exp])
			except lcompare.error as e:
				self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_call_2_bin(self):
		for nn in ((1,1,1), (5,5,5), (3,9,5,2,8)):
			totalrows = reduce(operator.mul, nn)
			exp = [reduce(operator.add, range(n))*(totalrows/n) for n in nn]
			cntcmd = ['cnt'] + [str(n) for n in nn]
			totcmd = ['total', '-of']
			res = self._runit([cntcmd, totcmd], ['run p1','run p2'])
			resl = struct.unpack('f'*len(exp), res)
			try: lcompare.lcompare(resl, exp)
			except lcompare.error as e:
					self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_call_many_text(self):
		data = b'''void polygon test 0 0 9  0 0 0  1 0 0  0 1 0'''
		exp = (['#', 'xform', '-s', 0.5],
				['#', 'xform', '-t', 0, 0, 2],
				['#', 'xform', '-ry', -20],
				['#', 'xform', '-rx', 30],
				['void', 'polygon', 'test'],
				[0], [0], [9],
				[0, 0, 1],
				[0.469846310393, 0, 1.17101007167],
				[-0.0855050358315, 0.433012701892, 1.23492315519])
		cmd0 = ['xform', '-rx', '30']
		cmd1 = ['xform', '-ry', '-20']
		cmd2 = ['xform', '-t', '0', '0', '2']
		cmd3 = ['xform', '-s', '0.5']
		cmdl = [cmd0, cmd1, cmd2, cmd3]
		resl = self._runit(cmdl, 'run many processes', indata=data,
				_in=PIPE, nl=True)
		try: lcompare.lcompare(resl, exp)
		except lcompare.error as e:
			print(resl, exp)
			self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_call_many_bin(self):
		exp = (15, 40, 45, 48, 0)
		cmd0 = ['cnt', '3', '4', '5', '1']
		cmd1 = ['histo', '1', '5']
		cmd2 = ['total', '-of']
		cmdl = [cmd0, cmd1, cmd2]
		res = self._runit(cmdl, 'run many processes')
		resl = struct.unpack('f'*len(exp), res)
		try: lcompare.lcompare(resl, exp)
		except lcompare.error as e:
			print(resl, exp)
			self.fail(('call_one_text -- ') +str(e))


# vi: set ts=4 sw=4 :
