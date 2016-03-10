

import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

_bindir = None

class CntTestCase(unittest.TestCase):
	def setUp(self):
		if _bindir:
			self.oldpath = os.environ['PATH']
			os.environ['PATH'] = _bindir

	def tearDown(self):
		if _bindir:
			os.environ['PATH'] = self.oldpath

	def test_1(self):
		cmd = 'cnt 5'
		res0 = os.popen(cmd).read()
		res = map(string.strip,string.split(res0, '\n'))
		exp = [0, 1, 2, 3, 4, '']
		try: lcompare.lcompare(res, exp)
		except lcompare.error, e: self.fail(str(e))

	def test_2(self):
		cmd = 'cnt 3 2'
		res0 = os.popen(cmd).read()
		res = map(string.split,string.split(res0, '\n'))
		exp = [[0,0], [0,1], [1,0], [1,1], [2,0], [2,1], []]
		try: lcompare.llcompare(res, exp)
		except lcompare.error, e: self.fail(str(e))

	def test_3(self):
		cmd = 'cnt 3 2 3'
		res0 = os.popen(cmd).read()
		res = map(string.split,string.split(res0, '\n'))
		exp = [[0,0,0],[0,0,1],[0,0,2],
			[0,1,0],[0,1,1],[0,1,2],
			[1,0,0],[1,0,1],[1,0,2],
			[1,1,0],[1,1,1],[1,1,2],
			[2,0,0],[2,0,1],[2,0,2],
			[2,1,0],[2,1,1],[2,1,2],
			[]]
		try: lcompare.llcompare(res, exp)
		except lcompare.error, e: self.fail(str(e))


def main(bindir=None):
	global _bindir
	_bindir=bindir
	support.run_case(CntTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
