
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class HistoTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath

	def test_histo(self):
		cmd = 'histo 1 6 5 < "%s"' % support.datafile('histo.dat')
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
			[1.5, 720, 240, 432, 1080, 270],
			[2.5, 720, 240, 432, 0, 270],
			[3.5, 0, 240, 432, 0, 270],
			[4.5, 0, 240, 432, 0, 270],
			[5.5, 0, 240, 0, 0, 270],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main():
	support.run_case(HistoTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
