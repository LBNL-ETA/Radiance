
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class XformTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath
		
	def test_xform(self):
		cmd = 'xform -e -s 3.14159 "%s"' % support.datafile('xform_1.dat')
		result = lcompare.split_rad(os.popen(cmd).read())
		expect = lcompare.split_radfile(support.datafile('xform_res1.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))


	def test_xform(self):
		cmd = 'xform -e -mx "%s"' % support.datafile('xform_2.dat')
		result = lcompare.split_rad(os.popen(cmd).read())
		expect = lcompare.split_radfile(support.datafile('xform_res2.dat'))
		try: lcompare.llcompare(result, expect, ignore_empty=1)
		except lcompare.error, e:
			self.fail('%s [%s]' % (str(e),cmd))

def main():
	support.run_case(XformTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
