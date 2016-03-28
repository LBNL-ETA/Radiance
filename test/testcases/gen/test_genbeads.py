# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import os
import unittest

from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class GenbeadsTestCase(unittest.TestCase, ProcMixin):

	def test_genbeads(self):
		cmdl = 'genbeads mymat myname 0 0 0 1 1 1 2 0 0 0 2 0 .1 .4'.split()
		try:
			proc = self.call_one(cmdl, 'call genbeads', out=PIPE,
				universal_newlines=True)
			raw = proc.stdout.read()
		except Error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmdl)))
		finally:
			proc.wait()
		result = lcompare.split_rad(raw)
		expect = [['mymat', 'sphere', 'myname.0'], [0], [0],
		[4, 0, 0, 0, 0.1],
		['mymat', 'sphere', 'myname.1'], [0], [0],
		[4, 0.36, 0.04, 0.104, 0.1],
		['mymat', 'sphere', 'myname.2'], [0], [0],
		[4, 0.651440715413, 0.167781092737, 0.365893348046, 0.1],
		['mymat', 'sphere', 'myname.3'], [0], [0],
		[4, 0.844350245496, 0.366600314978, 0.655866088042, 0.1],
		['mymat', 'sphere', 'myname.4'], [0], [0],
		[4, 0.960791445178, 0.643185551339, 0.897901825177, 0.1],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=True)
		except lcompare.error as e:
			self.fail('%s [%s]' % (str(e), self.qjoin(cmdl)))


# vi: set ts=4 sw=4 :
