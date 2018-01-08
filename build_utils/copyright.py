from __future__ import division, print_function, unicode_literals

import os
import sys
import builtins
if hasattr(builtins, 'raw_input'): # < Py3
	input = builtins.raw_input


def _get_ltextl():
	f = open(os.path.join('src','common','copyright.h'), 'r')
	ltl = f.readlines()
	f.close()
	ltl2 = []
	for line in ltl:
		line = line.strip()
		if line == '*/': line = ''
		elif line.find('$Id:') > -1: line = ''
		elif line and line[0] == '*':
			line = line[2:]
		elif line and line[1] == '*':
			line = line[3:]
		ltl2.append(line)
	return ltl2 + ['']

def _show_ltextl(ltextl, lines=23):
	llen = len(ltextl)
	for i in range(0, llen, lines):
		sys.stderr.write('\n'.join(ltextl[i:i+lines]))
		if i+lines < llen:
			input('\n[press <return> to continue] ')

def show_license():
	try:
		ltextl = _get_ltextl()
		_show_ltextl(ltextl)
		sys.stderr.write(
		'Do you understand and accept the terms of this agreement [n]?\n\n')
		answer = ''
		s =  'Please enter "yes" or "no", or use ^C to exit: '
		while answer not in ['y', 'ye', 'yes', 'n', 'no']:
			if answer: sys.stderr.write('invalid input "%s"\n' % answer)
			answer = input(s).lower()

		if answer[0] == 'y':
			return
		else:
			sys.stderr.write('\n*** Installation cancelled ***\n')
			sys.exit(1)
	except KeyboardInterrupt:
		sys.stderr.write('\n*** Installation cancelled ***\n')
		sys.exit(1)

# vi: set ts=4 sw=4 :

