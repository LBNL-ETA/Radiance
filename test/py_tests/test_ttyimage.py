

import os
import math
import string
import unittest

from unit_tools import support
from unit_tools import lcompare

class TtyimageTestCase(unittest.TestCase):
	def setUp(self):
		self.oldpath = os.environ['PATH']
		os.environ['PATH'] = os.path.abspath(support.BINDIR)

	def tearDown(self):
		os.environ['PATH'] = self.oldpath

	def test_ttyimage(self):
		'''We just do a few spot checks here'''
		picfile = support.datafile('Earth128.pic')
		cmd = 'ttyimage "%s"' % picfile
		res0 = os.popen(cmd).read()
		result = map(string.split,string.split(res0, '\n'))
		expect = [[0,
		['################################################################'
		 '################################################################']],
		[7,
		['#########################@%,,.?++&%%###$&###############@&:.....'
		 '.......,,.......,,.,,;..+?,...,.:*+.:&#########@################']],
		[23,
		['.......;,:.....,++*+?++++;+;:,::,..,,;+;;+...................,;,'
		 '..,:;+::+;;:;;:;;;;;:;+;;;;;;;;;:;+;;;;;%%$@%$&%#?.....,#%......']],
		[54,
		['................................,,,.......................*%?$@#'
		 '###########@$%%%;.........:?:.....,?+;:....;,...................']],
		[99,
		['.....................................+++........................'
		 '................................................................']],
		 ]

		for l in expect:
			self.assertEqual(result[l[0]], l[1],
			'%s : %s != %s [line %s]' % (cmd,result[l[0]],l[1], l[0]))

def main():
	support.run_case(TtyimageTestCase)

if __name__ == '__main__':
	main()

# vi: set ts=4 sw=4 :
