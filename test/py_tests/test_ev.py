
import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class EvTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath

	ltest = [
		['2.3 + 5 * 21.7 / 1.43', [2.3 + 5 * 21.7 / 1.43]],
		['if(1, 1, 0)', [1]],
		['if(0, 1, 0)', [0]],
		['select(3, 1, 2, 3, 4, 5)', [3]],
		['floor(5.743)', [5]],
		['ceil(5.743)', [6]],
		['sqrt(7.4)', [math.sqrt(7.4)]],
		['exp(3.4)', [math.exp(3.4)]],
		['log(2.4)', [math.log(2.4)]],
		['log10(5.4)', [math.log10(5.4)]],
		['sin(.51)', [math.sin(.51)]],
		['cos(.41)', [math.cos(.41)]],
		['tan(0.77)', [math.tan(0.77)]],
		['asin(0.83)', [math.asin(0.83)]],
		['acos(.94)', [math.acos(.94)]],
		['atan(0.22)', [math.atan(0.22)]],
		['atan2(0.72, 0.54)', [math.atan2(0.72, 0.54)]],
	]

	def test_singleres(self):
		for expr, expect in self.ltest:
			cmd = 'ev "%s"' % expr
			result = [string.strip(os.popen(cmd).read())]
			try: lcompare.lcompare(result, expect)
			except lcompare.error, e:
				self.fail('%s [%s]' % (str(e),cmd))

	def test_multipleres(self):
		pass # XXX implement

def main():
	support.run_case(EvTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
