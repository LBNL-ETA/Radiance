
import os
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
		cmd = 'histo -0.5 8.5 9 < "%s"' % support.datafile('histo.dat')
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
			[0, 720, 240, 432, 1080, 270],
			[1, 720, 240, 432, 1080, 270],
			[2, 720, 240, 432,    0, 270],
			[3,   0, 240, 432,    0, 270],
			[4,   0, 240, 432,    0, 270],
			[5,   0, 240,   0,    0, 270],
			[6,   0, 240,   0,    0, 270],
			[7,   0, 240,   0,    0, 270],
			[8,   0, 240,   0,    0,   0],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_histo_c(self):
		cmd = 'histo -c -0.5 8.5 9 < "%s"' % support.datafile('histo.dat')
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
			[-0.5,   0,    0,    0,    0,    0],
			[0.5,  720,  240,  432, 1080,  270],
			[1.5, 1440,  480,  864, 2160,  540],
			[2.5, 2160,  720, 1296, 2160,  810],
			[3.5, 2160,  960, 1728, 2160, 1080],
			[4.5, 2160, 1200, 2160, 2160, 1350],
			[5.5, 2160, 1440, 2160, 2160, 1620],
			[6.5, 2160, 1680, 2160, 2160, 1890],
			[7.5, 2160, 1920, 2160, 2160, 2160],
			[8.5, 2160, 2160, 2160, 2160, 2160],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_histo_p(self):
		cmd = 'histo -p -0.5 8.5 9 < "%s"' % support.datafile('histo.dat')
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
			[0, 33.333333, 11.111111, 20.0, 50.0, 12.5],
			[1, 33.333333, 11.111111, 20.0, 50.0, 12.5],
			[2, 33.333333, 11.111111, 20.0,  0.0, 12.5],
			[3,       0.0, 11.111111, 20.0,  0.0, 12.5],
			[4,       0.0, 11.111111, 20.0,  0.0, 12.5],
			[5,       0.0, 11.111111,  0.0,  0.0, 12.5],
			[6,       0.0, 11.111111,  0.0,  0.0, 12.5],
			[7,       0.0, 11.111111,  0.0,  0.0, 12.5],
			[8,       0.0, 11.111111,  0.0,  0.0,  0.0],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_histo_pc(self):
		cmd = 'histo -p -c -0.5 8.5 9 < "%s"' % support.datafile('histo.dat')
		raw = os.popen(cmd).read()
		result = map(string.split,string.split(raw,'\n'))
		expect = [
			[-0.5,      0.0,       0.0,   0.0,   0.0,   0.0],
			[0.5, 33.333333, 11.111111,  20.0,  50.0,  12.5],
			[1.5, 66.666667, 22.222222,  40.0, 100.0,  25.0],
			[2.5,     100.0, 33.333333,  60.0, 100.0,  37.5],
			[3.5,     100.0, 44.444444,  80.0, 100.0,  50.0],
			[4.5,     100.0, 55.555556, 100.0, 100.0,  62.5],
			[5.5,     100.0, 66.666667, 100.0, 100.0,  75.0],
			[6.5,     100.0, 77.777778, 100.0, 100.0,  87.5],
			[7.5,     100.0, 88.888889, 100.0, 100.0, 100.0],
			[8.5,     100.0,     100.0, 100.0, 100.0, 100.0],
		]
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main():
	support.run_case(HistoTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
