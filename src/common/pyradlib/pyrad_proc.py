# -*- coding: utf-8 -*-
''' pyrad_proc.py - Process and pipeline management for Python Radiance scripts
2016 - Georg Mischler

Use as:
	from pyradlib.pyrad_proc import PIPE, Error, ProcMixin

For a single-file installation, include the contents of this file
at the same place (minus the __future__ import below).
'''
from __future__ import division, print_function, unicode_literals

import sys
import subprocess
PIPE = subprocess.PIPE


class Error(Exception): pass


class ProcMixin():
	'''Process and pipeline management for Python Radiance scripts
	'''
	def raise_on_error(self, actstr, e):
		if hasattr(e, 'strerror'): eb = e.strerror
		elif isinstance(e, self._strtypes): eb = e
		else: eb = e
		if isinstance(eb, type(b'')):
			estr = eb.decode(encoding='utf-8', errors='ignore')
		else: estr = eb
		raise Error('Unable to %s - %s' % (actstr, estr)) #

	def __configure_subprocess(self):
		'''Prevent subprocess module failure in frozen scripts on Windows.
		   Prevent console windows from popping up when not console based.
		   Make sure we use the version-specific string types.
		'''
		# On Windows, sys.stdxxx may not be available when:
		# - built as *.exe with "pyinstaller --noconsole"
		# - invoked via CreateProcess() and stream not redirected
		try:
			sys.__stdin__.fileno()
			self._stdin = sys.stdin
		except: self._stdin = PIPE
		try:
			sys.__stdout__.fileno()
			self._stdout = sys.stdout
		except: self._stdout = PIPE
		try:
			sys.__stderr__.fileno()
			self._stderr = sys.stderr
		# keep subprocesses from opening their own console.
		except: self._stderr = PIPE
		if hasattr(subprocess, 'STARTUPINFO'):
			si = subprocess.STARTUPINFO()
			si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
			self._pipeargs = {'startupinfo':si}
		else: self._pipeargs = {}
		# type names vary between Py2.7 and 3.x
		self._strtypes = (type(b''), type(u''))
		# private attribute to indicate established configuration
		self.__proc_mixin_setup = True

	def qjoin(self, sl):
		'''Join a list with quotes around each element containing whitespace.
		We only use this to display command lines on sys.stderr, the actual
		Popen() calls are made with the original list.
		'''
		def _q(s):
			if ' ' in s or '\t' in s or ';' in s:
				return "'" + s + "'"
			return s
		return  ' '.join([_q(s) for s in sl])

	def __parse_args(self, _in, out):
		try: self.__proc_mixin_setup
		except AttributeError: self.__configure_subprocess()
		instr = ''
		if _in == PIPE:
			stdin = _in
		elif isinstance(_in, self._strtypes):
			if self.donothing: stdin = None
			else: stdin = open(_in, 'rb')
			instr = ' < "%s"' % _in
		elif hasattr(_in, 'read'):
			stdin = _in
			instr = ' < "%s"' % _in.name
		else: stdin = self._stdin
		outstr = ''
		if out == PIPE:
			stdout = out
		elif isinstance(out, self._strtypes):
			if self.donothing: stdout = None
			else: stdout = open(out, 'wb')
			outstr = ' > "%s"' % out
		elif hasattr(out, 'write'):
			stdout = out
			outstr = ' > "%s"' % out.name
		else: stdout = self._stdout
		return stdin, stdout, instr, outstr

	def call_one(self, cmdl, actstr, _in=None, out=None,
			universal_newlines=False):
		'''Create a single subprocess, possibly with an incoming and outgoing
		pipe at each end.
		- cmdl
		  A list of strings, leading with the name of the executable (without
		  the *.exe suffix), followed by individual arguments.
		- actstr
		  A text string of the form "do something".
		  Used in verbose mode as "### do someting\\n### [command line]"
		  Used in error messages as "Scriptname: Unable to do something".
		- _in / out
		  What to do with the input and output pipes of the process:
		  * a filename as string
		    Open file and use for reading/writing.
		  * a file like object
		    Use for reading/writing
		  * PIPE
		    Pipe will be available in returned object for reading/writing.
		  * None (default)
		    System stdin/stdout will be used if available
		If _in or out is a PIPE, the caller should call p.wait() on the
		returned Popen instance after writing to and closing it.
		'''
		stdin, stdout, instr, outstr = self.__parse_args(_in, out)
		if getattr(self, 'verbose', None):
			sys.stderr.write('### %s \n' % actstr)
			sys.stderr.write(self.qjoin(cmdl) + instr + outstr + '\n')
		if not getattr(self, 'donothing', None):
			try: p = subprocess.Popen(cmdl,
					stdin=stdin, stdout=stdout, stderr=self._stderr,
					universal_newlines=universal_newlines, **self._pipeargs)
			except Exception as e:
				self.raise_on_error(actstr, e)
			if stdin != PIPE and stdout != PIPE:
				# caller needs to wait after reading or writing (else deadlock)
				res = p.wait()
				if res != 0:
					self.raise_on_error(actstr,
							'Nonzero exit (%d) from command [%s].'
							% (res, self.qjoin(cmdl)+instr+outstr+'\n'))
			return p

	def call_two(self, cmdl_1, cmdl_2, actstr_1, actstr_2, _in=None, out=None,
			universal_newlines=False):
		'''Create two processes, chained via a pipe, possibly with an incoming
		and outgoing pipe at each end.
		Returns a tuple of two Popen instances.
		Arguments are equivalent to call_one(), with _in and out applying
		to the ends of the chain.
		If _in or out is PIPE, the caller should call p.wait() on both
		returned popen instances after writing to and closing the first on .
		'''
		stdin, stdout, instr, outstr = self.__parse_args(_in, out)
		if getattr(self, 'verbose', None):
			sys.stderr.write('### %s \n' % actstr_1)
			sys.stderr.write('### %s \n' % actstr_2)
			sys.stderr.write(self.qjoin(cmdl_1) + instr + ' | ')
		if not getattr(self, 'donothing', None):
			try: p1 = subprocess.Popen(cmdl_1,
					stdin=stdin, stdout=PIPE, stderr=self._stderr,
					**self._pipeargs)
			except Exception as e:
				self.raise_on_error(actstr_1, e)
		if getattr(self, 'verbose', None):
			sys.stderr.write(self.qjoin(cmdl_2) + outstr + '\n')
		if not getattr(self, 'donothing', None):
			try:
				p2 = subprocess.Popen(cmdl_2,
						stdin=p1.stdout, stdout=stdout, stderr=self._stderr, 
						universal_newlines=universal_newlines, **self._pipeargs)
				p1.stdout.close()
			except Exception as e:
				self.raise_on_error(actstr_2, e)
			if stdin != PIPE and stdout != PIPE:
				# caller needs to wait after reading or writing (else deadlock)
				res = p1.wait()
				if res != 0:
					self.raise_on_error(actstr_1,
							'Nonzero exit (%d) from command [%s].'
							% (res, self.qjoin(cmdl_1)))
				res = p2.wait()
				if res != 0:
					self.raise_on_error(actstr_2,
							'Nonzero exit (%d) from command [%s].'
							% (res, self.qjoin(cmdl_2)))
			return p1, p2

	def call_many(self, cmdlines, actstr, _in=None, out=None,
			universal_newlines=False):
		'''Create a series of N processes, chained via pipes, possibly with an
		incoming and outgoing pipe at each end.
		Returns a tuple of N subprocess.Popen instances.
		Depending on the values of _in and out, the first and last may be
		available to write to or read from respectively.
		Most arguments are equivalent to call_one(), with
		- cmdlines
		  a list of N command argument lists
		- _in / out
		  applying to the ends of the chain.
		If _in or out is PIPE, the caller should call p.wait() on all
		returned Popen instances after writing to and closing the first one.
		'''
		if len(cmdlines) == 1:
			# other than direct call_one(), this returns a one-item tuple!
			return (self.call_one(cmdlines[0], actstr, _in=_in, out=out,
					universal_newlines=universal_newlines),)
		stdin, stdout, instr, outstr = self.__parse_args(_in, out)
		procs = []
		if getattr(self, 'verbose', None):
			sys.stderr.write('### %s \n' % actstr)
			sys.stderr.write(self.qjoin(cmdlines[0]) + instr + ' | ')
		if not getattr(self, 'donothing', None):
			try:
				prevproc = subprocess.Popen(cmdlines[0], stdin=stdin,
						stdout=PIPE, stderr=self._stderr, **self._pipeargs)
				procs.append(prevproc)
			except Exception as e:
				self.raise_on_error(actstr, e)

		for cmdl in cmdlines[1:-1]:
			if getattr(self, 'verbose', None):
				sys.stderr.write(self.qjoin(cmdl) + ' | ')
			if not getattr(self, 'donothing', None):
				try:
					nextproc = subprocess.Popen(cmdl, stdin=prevproc.stdout,
							stdout=PIPE, stderr=self._stderr, **self._pipeargs)
					procs.append(nextproc)
					prevproc.stdout.close()
					prevproc = nextproc
				except Exception as e:
					self.raise_on_error(actstr, e)

		if getattr(self, 'verbose', None):
			sys.stderr.write(self.qjoin(cmdlines[-1]) + outstr + '\n')
		if not getattr(self, 'donothing', None):
			try:
				lastproc = subprocess.Popen(cmdlines[-1], stdin=prevproc.stdout,
						stdout=stdout, stderr=self._stderr,
						universal_newlines=universal_newlines, **self._pipeargs)
				prevproc.stdout.close()
				procs.append(lastproc)
				prevproc.stdout.close()
			except Exception as e:
				self.raise_on_error(actstr, e)

			if stdin != PIPE and stdout!= PIPE:
				# caller needs to wait after reading or writing (else deadlock)
				for proc, cmdl in zip(procs, cmdlines):
					res = proc.wait()
					if res != 0:
						self.raise_on_error(actstr,
								'Nonzero exit (%d) from command [%s].'
								% (res, self.qjoin(cmdl)))
			return procs


### end of proc_mixin.py
