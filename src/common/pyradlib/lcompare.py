# -*- coding: utf-8 -*-
''' Text comparison functions for Radiance unit testing.

The comparison allows for differences in whitespace, so we need to
split the text into tokens first. Tokens are then converted into an
appropriate data type, so that we can do an "almost equal" comparison
for floating point values, ignoring binary rounding errors.
'''
from __future__ import division, print_function, unicode_literals
__all__ = ['error', 'lcompare', 'llcompare',
		'split_headers', 'split_rad', 'split_radfile']

# Py2.7/3.x compatibility
try: from itertools import (izip_longest as zip_longest, chain,
		ifilter as filter, izip as zip)
except ImportError:
	from itertools import zip_longest, chain

class error(Exception): pass

# internal functions
def _icompare(itest, iref):
	'''compare ints (not public)'''
	if type(itest) == str:
		iftest = int(itest)
	else: iftest = itest
	if iftest == iref: return 1
	return 0

def _fcompare(ftest, fref):
	'''compare floats (not public)'''
	FUZZ = 0.0000001 # XXX heuristically determined (quite low?)
	if type(ftest) == str:
		fftest = float(ftest)
	else: fftest = ftest
	if (fftest < (fref + FUZZ)) and (fftest > (fref - FUZZ)):
		return True
	return False

def _typify_token(t):
	'''return the token as int resp. float if possible (not public)'''
	try: return int(t)
	except ValueError: pass
	try: return float(t)
	except ValueError: pass
	return t


# public comparison functions

def lcompare(ltest, lref):
	'''Compare a list/iterator of tokens.
	   Raise an error if there are intolerable differences.
	   The reference tokens in lref should already be of the correct type.
	'''
	i = 0
	for ttest, tref in zip_longest(ltest, lref, fillvalue=False):
		if ttest is False:
			raise error('List comparision failed: Fewer tokens than expected'
				' (%d  !=  >= %d)' % (i, i+1))
		if tref is False:
			raise error('List comparision failed: More tokens than expected'
				' (>= %d  !=  %d)' % (i+1, i))
		if type(tref) == str and tref != ttest:
			print(tref, ttest)
			raise error('String token comparison failed: "%s" != "%s"'
					% (ttest, tref))
		elif type(tref) == int and not _icompare(ttest, tref):
			print((ttest, tref))
			raise error('Int token comparison failed: %s != %s' % (ttest, tref))
		elif type(tref) == float and not _fcompare(ttest, tref):
			raise error('Float token comparison failed: %s != %s'
					% (ttest, tref))
		i += 1

def llcompare(lltest, llref, ignore_empty=False, _recurse=[]):
	'''Compare a list/iterator of lists/iterators of tokens recursively.
	   Raise an error if there are intolerable differences.
	   The reference tokens in lref should already be of the correct type.
	   If ignore_empty is true, empty lines are not included in the comparison.
	   The _recurse argument is only used internally.
	'''
	if ignore_empty:
		lltest = filter(None, lltest)
		llref = filter(None, llref)
	i = 0
	for ltest, lref in zip_longest(lltest, llref, fillvalue=False):
		i += 1
		if ltest is False:
			raise error('List comparision failed: Fewer entries than expected'
				' (%d  !=  >= %d)' % (i, i+1))
		if lref is False:
			raise error('List comparision failed: More entries than expected'
				' (>= %d  !=  %d)' % (i+1, i))
		if lref and not isinstance(lref, str):
			if hasattr(lref, '__getitem__'):
				rfirst = lref[0]
			elif hasattr(lref, 'next') or hasattr(lref, '__next__'):
				rfirst = next(lref) # "peek" at first
				lref = chain([rfirst], lref) # "push" back
			else: rfirst = None
			if isinstance(rfirst, str):
				rfirst = None
			if hasattr(rfirst, '__iter__') or isinstance(rfirst, (list, tuple)):
				return llcompare(ltest, lref,
						_recurse=_recurse + [i], ignore_empty=ignore_empty)
		try: lcompare(ltest, lref)
		except TypeError:
			print(ltest, lref)
			raise
		except error as e:
			if _recurse:
				raise error('%s (line %s)' % (str(e), _recurse + [i + 1]))
			else: raise error('%s (line %d)' % (str(e), i + 1))
	
def split_headers(s):
	'''split Radiance file headers (eg. the output of getinfo).
	   Return a list of lists of tokens suitable for llcompare().'''
	ll = [ss.strip() for ss in s.split('\n')]
	nll = []
	for l in ll:
		parts = l.split('=', 1)
		if len(parts) == 2:
			left = [_typify_token(s) for s in parts[0].split()]
			right = [_typify_token(s) for s in parts[1].split()]
			nll.append(left + ['='] + right)
		else: nll.append([_typify_token(s) for s in l.split()])
	return nll

def split_rad(s):
	'''Split the contents of a scene description string.
	   Return a list of list of tokens suitable for llcompare().'''
	ll = [ss.strip() for ss in s.split('\n')]
	return [[_typify_token(s) for s in l.split()] for l in ll]

def split_radfile(fn):
	'''Split the contents of a file containing a scene description.
	   Return a list of list of tokens suitable for llcompare().'''
	with open(fn, 'r') as f:
		return split_rad(f.read())
	
