''' Text comparison functions for Radiance unit testing.

This allows differences in whitespace, which is why the text
corpora are split into tokens first.
Tokens are then converted into an appropriate data type, so
that floating point items will still be considered correct
even if they are slightly different, eg. as a consequence of
binary rounding errors.
'''

import string
import types

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
	FUZZ = 0.0000001 # XXX heuristically determined
	if type(ftest) == str:
		fftest = float(ftest)
	else: fftest = ftest
	if (fftest < (fref + FUZZ)) and (fftest > (fref - FUZZ)):
		return 1
	return 0

def _typify_token(t):
	'''return the token as int resp. float if possible (not public)'''
	try: return int(t)
	except ValueError: pass
	try: return float(t)
	except ValueError: pass
	return t


# public comparison functions

def lcompare(ltest, lref):
	'''compare a list of tokens
		raise an error if there are intolerable differences
		the reference tokens in lref should already be of the correct type.
	'''
	if len(ltest) != len(lref):
		raise error, ('List comparision failed: Different number of tokens'
			' (%d, %d)' % (len(ltest), len(lref)))
	for i in range(len(lref)):
		tref = lref[i]
		ttest = ltest[i]
		if type(tref) == str and tref != ttest:
			raise error, 'Token comparison failed: "%s" != "%s"' % (ttest, tref)
		elif type(tref) == int and not _icompare(ttest, tref):
			raise error, 'Token comparison failed: %s != %s' % (ttest, tref)
		elif type(tref) == float and not _fcompare(ttest, tref):
			raise error, 'Token comparison failed: %s != %s' % (ttest, tref)

def llcompare(lltest, llref, ignore_empty=0, recurse=[]):
	'''compare a list of lists of tokens recursively
		raise an error if there are intolerable differences
		the reference tokens in lref should already be of the correct type.
		if ignore_empty is true, empty lines are not included in the comparison
		the recurse argument is only used internally
	'''
	if ignore_empty:
		lltest = filter(None, lltest)
		llref = filter(None, llref)
	if len(lltest) != len(llref):
		raise error, 'Comparision failed: Different number of lines (%d,%d)' %(
			len(lltest), len(llref))
	for i in range(len(llref)):
		if llref[i]:
			rtype = type(llref[i][0])
			if rtype == list or rtype == tuple:
				return llcompare(lltest[i], llref[i],
						recurse=recurse.append(i), ignore_empty=ignore_empty)
		try: lcompare(lltest[i], llref[i])
		except error, e:
			if recurse:
				raise error, '%s (line %s)' % (str(e), recurse.append(i + 1))
			else: raise error, '%s (line %d)' % (str(e), i + 1)
	
def split_headers(s):
	'''split Radiance file headers
		return a list of lists of tokens suitable for llcompare()
		this is useful to check the output of getinfo'''
	ll = map(string.strip,string.split(s, '\n'))
	nll = []
	for l in ll:
		parts = string.split(l, '=', 1)
		if len(parts) == 2:
			left = map(_typify_token, string.split(parts[0]))
			right = map(_typify_token, string.split(parts[1]))
			nll.append(left + ['='] + right)
		else: nll.append(map(_typify_token, string.split(l)))
	return nll

def split_rad(s):
	'''Split the contents of a scene description string
		return a list of list of tokens suitable for llcompare()'''
	ll = map(string.strip,string.split(s, '\n'))
	nll = []
	for l in ll:
		nll.append(map(_typify_token, string.split(l)))
	return nll

def split_radfile(fn):
	'''Split the contents of a file object containing a scene description
		return a list of list of tokens suitable for llcompare()'''
	f = open(fn, 'r')
	ll = split_rad(f.read())
	f.close()
	return ll
	
