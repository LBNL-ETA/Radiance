# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import unittest

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class TtyimageTestCase(unittest.TestCase, ProcMixin):

	def _runit(self, cmd):
		try:
			proc = self.call_one(cmd, 'call ttyimage', out=PIPE,
					universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmd)))
		finally:
			proc.wait()
		return [s.split() for s in raw.split('\n')]

	def test_ttyimage(self):
		# We just do a few spot checks here
		picfile = ts.datafile('Earth128.pic')
		cmd = ['ttyimage', picfile]
		result = self._runit(cmd)
		expect = [[0,
		['################################################################'
		 '################################################################']],
		[7,
		['#########################@%,,.?++&%%###$&###############@&:.....'
		 '.......,,.......,,.,,;..+?,...,.:*+.:&#########@################']],
		[23,
		['.......;,:.....,++*+?++++;+;:,::,..,,;+;;+...................,;,'
		 '..,:;+::+;;:;;:;;;;;:;+;;;;;;;;;:;+;;;;;%%$@%$&%#?.....,#%......']],
		[54,
		['................................,,,.......................*%?$@#'
		 '###########@$%%%;.........:?:.....,?+;:....;,...................']],
		[99,
		['.....................................+++........................'
		 '................................................................']],
		 ]

		for l in expect:
			self.assertEqual(result[l[0]], l[1],
			'%s : %s != %s [line %s]' % (cmd,result[l[0]],l[1], l[0]))


# vi: set ts=4 sw=4 :
