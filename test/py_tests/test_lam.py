
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class LamTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath

	def test_lam(self):
		dat_de = support.datafile('lam_de.dat')
		dat_en = support.datafile('lam_en.dat')
		cmd = 'rlam -t: "%s" "%s"' % (dat_de, dat_en)
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
		['eins:one'],
		['zwei:two'],
		['drei:three'],
		['vier:four'],
		['fuenf:five'],
	]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main():
	support.run_case(LamTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
