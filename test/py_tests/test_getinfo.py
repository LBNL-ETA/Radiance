

import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class GetinfoTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath


	def test_getinfo(self):
		picfile = support.datafile('Earth128.pic')
		cmd = 'getinfo "%s"' % picfile
		res0 = os.popen(cmd).read()
		result = lcompare.split_headers(res0)
		exps = '''%s:
        Xim format conversion by:
        FORMAT=32-bit_rle_rgbe
        pfilt -e 2 -x 512 -y 512 -p 1 -r .67
        EXPOSURE=4.926198e+00
        normpat
        pfilt -1 -e .2
        EXPOSURE=2.000000e-01
        pfilt -x 128 -y 128
        PIXASPECT=0.500000
        EXPOSURE=2.571646e+00''' % picfile
		expect = lcompare.split_headers(exps)
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

	def test_getinfo_d(self):
		picfile = support.datafile('Earth128.pic')
		cmd = 'getinfo -d "%s"' % picfile
		res0 = os.popen(cmd).read()
		result = map(string.split,string.split(res0, '\n'))
		exps = '''%s: -Y 128 +X 128''' % picfile
		expect = map(string.split, string.split(exps, '\n'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main():
	support.run_case(GetinfoTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
