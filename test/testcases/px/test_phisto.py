# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import os
import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class PhistoTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd, _in=None):
		failmsg = None
		proc = None
		try:
			proc = self.call_one(cmd, 'test phisto',
					_in=_in, out=PIPE, universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			#self.fail(failmsg)
			failmsg = str(e)
		finally:
			if proc:
				res = proc.wait()
				if res != 0:
					failmsg = 'Non-zero (%d) exit from %s' % (res, str(cmd))
		if failmsg:
			self.fail(failmsg)
		return lcompare.split_rad(raw)

	def test_phisto(self):
		if os.name == 'nt': phisto = 'phisto'
		else: phisto = 'phisto'
		hgradpic = ts.datafile('gradients', 'h_gradient.hdr')
		vgradpic = ts.datafile('gradients', 'v_gradient.hdr')
		datafn = ts.datafile('gradients', 'gradient.histo')
		with open(datafn, 'r') as df:
			dtxt = df.read()
		expect = lcompare.split_rad(dtxt)
		for picfn in (hgradpic, vgradpic):
			hcmd = [phisto, picfn]
			err_msg = None
			try:
				result = self._runit(hcmd)
				lcompare.llcompare(result, expect)
			except AssertionError as e:
				with open(picfn, 'rb') as picf:
					hcmd = [phisto]
					result = self._runit(hcmd, _in=picf)
					try:
						lcompare.llcompare(result, expect)
						err_msg = 'Phisto fails with spaces in file name paths.'
					except Exception as e:
						err_msg = str(e)
			except Exception as e:
				err_msg = str(e)
			if err_msg:
				self.fail(err_msg)



# vi: set ts=4 sw=4 :
