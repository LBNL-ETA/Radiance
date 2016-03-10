
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

_bindir = None

class XformTestCase(unittest.TestCase):
	def setUp(self):
		if _bindir:
			self.oldpath = os.environ['PATH']
			os.environ['PATH'] = _bindir

	def tearDown(self):
		if _bindir:
			os.environ['PATH'] = self.oldpath
		
	def test_xform_e_s(self):
		cmd = 'xform -e -s 3.14159 "%s"' % support.datafile('xform_1.dat')
		result = lcompare.split_rad(os.popen(cmd).read())
		expect = lcompare.split_radfile(support.datafile('xform_res1.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))


	def test_xform_e_mx(self):
		cmd = 'xform -e -mx "%s"' % support.datafile('xform_2.dat')
		result = lcompare.split_rad(os.popen(cmd).read())
		expect = lcompare.split_radfile(support.datafile('xform_res2.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main(bindir=None):
	global _bindir
	_bindir=bindir
	support.run_case(XformTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
