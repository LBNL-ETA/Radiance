
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class GenbeadsTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath

	def test_genbeads(self):
		cmd = 'genbeads mymat myname 0 0 0 1 1 1 2 0 0 0 2 0 .1 .4'
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
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
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))


def main():
	support.run_case(GenbeadsTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
