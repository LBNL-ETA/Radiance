# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals

import struct
import operator
import unittest
from functools import reduce

import testsupport as ts
from pyradlib import lcompare
from pyradlib.pyrad_proc import PIPE, Error, ProcMixin


class LcompareTestCase(unittest.TestCase, ProcMixin):

	def test_lc_llcompare(self):
		# with values higher than 44721, total will return a float in e-format.
		for xc, data, exp in (
				(  None,
				   ( ('abcde', 'fgh', '1234', '56.789'),
					 ('xyz', '432', '987.65432')
				   ),
				   ( ('abcde', 'fgh', 1234, 56.789),
					 ('xyz', 432, 987.65432)
				   ),
				),
				(  None,
				   ( (('a', 'b', 'c'),('d', 'e', 'f')),
					 (('0', '1', '2', '3'),('1.1', '2.2', '3.000e-03')),
				   ),
				   ( (('a', 'b', 'c'),('d', 'e', 'f')),
					 (range(4),(1.1, 2.2, 0.003)),
				   ),
				),
				(  lcompare.error, # top level length
				   (('a', 'b', 'c'),('c','d'),('e','f'),),
				   (('a', 'b', 'c'),('c','d'),),
				),
				(  lcompare.error, # top level length
				   (('a', 'b', 'c'),('c','d'),),
				   (('a', 'b', 'c'),('c','d'),('e','f'),),
				),
				(  lcompare.error, # 2nd level length
				   (('a', 'b', 'c'),('c','d'),('e','f','g'),),
				   (('a', 'b', 'c'),('c','d'),('e','f'),),
				),
				(  lcompare.error, # 2nd level length
				   (('a', 'b', 'c'),('c','d'),('e','f'),),
				   (('a', 'b', 'c'),('c','d'),('e','f','g'),),
				),
				(  lcompare.error, # string diff
				   (('a', 'b', 'c'),('c','d'),('e','f','g'),),
				   (('a', 'b', 'c'),('c','d'),('e','f','h'),),
				),
				(  lcompare.error, # int diff
				   (('a', 'b', 'c'),('c','d'),('1','2','3'),),
				   (('a', 'b', 'c'),('c','d'),( 1,  2,  4),),
				),
				(  lcompare.error, # float diff
				   (('a', 'b', 'c'),('c','d'),('1.1','2.2','3.3'),),
				   (('a', 'b', 'c'),('c','d'),( 1.1,  2.2,  3.4),),
				),
				(  lcompare.error, # exponent diff
				   (('a', 'b', 'c'),('c','d'),('1.1','2.2','3.0000e-02'),),
				   (('a', 'b', 'c'),('c','d'),( 1.1,  2.2,  0.003),),
				),
				(  lcompare.error, # fuzzy compare
				   (('a', 'b', 'c'),('c','d'),('1.1','2.2','3.00000003'),),
				   (('a', 'b', 'c'),('c','d'),( 1.1,  2.2,  3.0000003),),
				),
				(  None, # fuzzy compare 
				   (('a', 'b', 'c'),('c','d'),('1.1','2.2','3.000000003'),),
				   (('a', 'b', 'c'),('c','d'),( 1.1,  2.2,  3.00000003),),
				),
				):
			if xc:
				self.assertRaises(xc, lcompare.llcompare, data, exp)
			else:
				try: lcompare.llcompare(data, exp)
				except lcompare.error as e:
					self.fail(('call_one_text ') +str(e))

	def test_lc_split_headers(self):
		htxt = '''example.hdr:
		Xim format conversion by:
		FORMAT=32-bit_rle_rgbe
		pfilt -e 2 -x 512 -y 512 -p 1 -r .67
		EXPOSURE=4.926198e+00
		normpat
		pfilt -1 -e .2
		EXPOSURE=2.000000e-01
		pfilt -x 128 -y 128
		PIXASPECT=0.500000
		EXPOSURE=2.571646e+00'''
		res = lcompare.split_headers(htxt)
		exp = (
				('', ('example.hdr:',)),
				('\t\t', ('Xim','format','conversion','by:')),
				('\t\t', 'FORMAT', '=', ('32-bit_rle_rgbe',)),
				('\t\t',
				('pfilt','-e','2','-x','512','-y','512','-p','1','-r','.67')),
				('\t\t', 'EXPOSURE', '=', ('4.926198e+00',)),
				('\t\t', ('normpat',)),
				('\t\t', ('pfilt','-1','-e','.2',)),
				('\t\t', 'EXPOSURE', '=', ('2.000000e-01',)),
				('\t\t', ('pfilt','-x','128','-y','128',)),
				('\t\t', 'PIXASPECT', '=', ('0.500000',)),
				('\t\t', 'EXPOSURE', '=', ('2.571646e+00',)),
			  )
		try: lcompare.llcompare(res, exp)
		except lcompare.error as e:
				self.fail(('call_one_text ') +str(e))

	def test_lc_split_radfile(self):
		df = ts.datafile('window_src.rad')
		exp = ([['#'],
				['#', 'A', 'plain', 'old', 'glass', 'window'],
				['#'],
				[], ['void', 'light', 'window_light'],
				[0], [0], [3, 1, 1, 1],
				[], ['window_light', 'polygon', 'window'],
				[0], [0], [12], [23.5, 43, 30], [23.5, 26, 30],
				[-23.5, 26, 30], [-23.5, 43, 30], []])
		resl = lcompare.split_radfile(df)
		try: lcompare.lcompare(resl, exp)
		except lcompare.error as e:
			print(resl, exp)
			self.fail(('call_one_text n=%d -- ' % n) +str(e))

	def test_lc_split_rad(self):
		df = ts.datafile('window_src.rad')
		exp = ([['#'],
				['#', 'A', 'plain', 'old', 'glass', 'window'],
				['#'],
				[], ['void', 'light', 'window_light'],
				[0], [0], [3, 1, 1, 1],
				[], ['window_light', 'polygon', 'window'],
				[0], [0], [12], [23.5, 43, 30], [23.5, 26, 30],
				[-23.5, 26, 30], [-23.5, 43, 30], []])
		with open(df) as f:
			res = f.read()
		resl = lcompare.split_rad(res)
		try: lcompare.lcompare(resl, exp)
		except lcompare.error as e:
			print(resl, exp)
			self.fail(('call_one_text n=%d -- ' % n) +str(e))


# vi: set ts=4 sw=4 :
