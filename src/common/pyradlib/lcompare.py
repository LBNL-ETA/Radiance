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

import re
import shlex
# Py2.7/3.x compatibility
try: from itertools import (izip_longest as zip_longest, chain,
		ifilter as filter, izip as zip)
except ImportError:
	from itertools import zip_longest, chain

class error(Exception): pass

_strtypes = (type(b''), type(u''))

# internal functions
def _icompare(itest, iref):
	'''compare ints (not public)'''
	if isinstance(itest, _strtypes):
		iftest = int(itest)
	else: iftest = itest
	if iftest == iref: return True
	return False

def _fcompare(ftest, fref):
	'''compare floats (not public)'''
	FUZZ = 0.0000001 # XXX heuristically determined (quite low?)
	if isinstance(ftest, _strtypes):
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
		if isinstance(tref, _strtypes) and tref != ttest:
			raise error('String token comparison failed: "%s" != "%s"'
					% (ttest, tref))
		elif type(tref) == int and not _icompare(ttest, tref):
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
		if lref and not isinstance(lref, _strtypes):
		
			if hasattr(lref, '__getitem__'):
				rfirst = lref[0]
			elif hasattr(lref, 'next') or hasattr(lref, '__next__'):
				rfirst = next(lref) # "peek" at first
				lref = chain([rfirst], lref) # "push" back
			else: rfirst = None
			if isinstance(rfirst, _strtypes):
				rfirst = None
			if hasattr(rfirst, '__iter__') or isinstance(rfirst, (list, tuple)):
				llcompare(ltest, lref,
						_recurse=_recurse + [i], ignore_empty=ignore_empty)
		try: lcompare(ltest, lref)
#		except TypeError:
#			print(ltest, lref)
#			raise
		except error as e:
			if _recurse:
				raise error('%s (line %s)' % (str(e), _recurse + [i + 1]))
			else: raise error('%s (line %d)' % (str(e), i + 1))
	
_HLPATS  = '(\s*)(?:([^=\s]*)\s*=\s*(.*)\s*|(.*)\s*)'
_hlpat = re.compile(_HLPATS)
def split_headers(s):
	'''split Radiance file headers (eg. the output of getinfo).
	   Return a list of lists of tokens suitable for llcompare().'''
	ll = s.split('\n')
	nll = []
	for l in ll:
		m = _hlpat.match(l)
		groups = m.groups()
		indent = groups[0]
		if groups[1]:
			left = groups[1]
			right = [_typify_token(s) for s in shlex.split(groups[2])]
			nll.append([indent, left, '='] + [right])
		else:
			full = [_typify_token(s) for s in shlex.split(groups[3])]
			nll.append([indent] + [full])
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
	
