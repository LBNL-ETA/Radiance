
Radiance Testing Framework
--------------------------

A toolkit to test all (eventually) components of the Radiance
synthetic image generation system for conformance to their
specification.


Limitations

For the moment, we use PyUnit to run our tests. This means that
we're restricted to test only complete programs, and not actual
units (since PyUnit was designed to test Python units, not C).
A C-level testing framework may be added later.

There's no good way to automatically test GUI programs like
rview. We have to rely on good human testers to check whether
those work correctly or not.


Requirements

You need a working installation of Python 2.1 (or newer) on your
system. The reason for this is that the PyUnit framework isn't
included with earlier versions. If you prefer to use an older
Python (back to 1.5.2), you can get PyUnit here, and install it
somewhere on your PYTHONPATH:
  http://pyunit.sourceforge.net/

Our testing framework currently assumes that the Radiance files
reside in the following local file tree (seen from the "test/"
subdirectory where this file resides):

  executables:    ../bin/*[.exe]
  support files:  ../lib/*
  data files:     ./test data/*

This is the location where the SCons build system places
everything, which means we're testing the software after building
but before installing. 
The space character in the name of the test data directory is
deliberate, because it is a design requirement that all our
executables can handle path names with spaces.


How to run tests

The simplest way to run tests is to use the SCons build system.
The file ray/INSTALL.scons explains the requirements and details.
Once you have SCons working, go to the ray directory and type

 $> scons test

which will automatically execute all available tests in the
correct environment.

You can also run the tests manually:

On unix systems, just type "run_all.py" in this directory to
test everything. If that file doesn't have execute rights, you
can supply it to the python interpreter as its single argument:
"python run_all.py". You can also run individual test suites from
the "py_tests" directory directly: "python test_getinfo.py".

On Windows, this should usually work as well. As an alternative,
use the "winrun.bat" script. WARNING: You need to change the
paths in this script to match your Python installation first.


What gets tested

There are several test groups, each containing a number of test
suites, each containing one or more tests. When running tests,
the name of the test groups and test suites will get printed to
the console, the latter with an "ok" if they succeeded.

If any test fails, there will be diagnostic output about the
nature of the failure, but the remaining tests will continue to
be executed. Note that several utility programs may be used to
access the results of other calculations, so if eg. getinfo is
broken, that may cause a number of seemingly unrelated tests to
fail as well.


How to report failures

If any of the tests fail on your platform, please report your
results (and as much ancilliary information about your system and
Radiance version as possible) to the radiance code development
mailing list on http://www.radiance-online.org/
The developers will then either try to fix the bug, or instruct
you on how to refine your testing to get more information about
what went wrong.


How to contribute test cases

The list of tests run is still very much incomplete, but will
hopefully grow quickly. You can contribute by creating tests too!
Please ask on the code development mailing list first, so that we
can avoid overlaps between the work of different contributors.

There are two classes of tests to be considered:

- Testing individual executables
  This means that an individual program like ev, xfom, or getinfo
  is tested with typical input data, and the output is compared
  against the expected result.

- Testing specific calculations
  This will mainly affect the actual simulation programs rpict
  and rtrace. For example, there should be a test suite for every
  material (and modifier) type, which uses rtrace to shoot a
  series of rays against a surface under varying angles, in order
  to verify material behaviour under different parameters. Tests
  of this kind may require a custom script.

Contributed tests can be of two kinds. In the simplest case, you
can contribute a small(!) set of test data, the command line(s)
used to run your tests on them, and a list of expected results.
Result comparisons are typically done in text form (by line).
If the result is a picture, we'll use ttyimage to pick out a few
scan lines for comparison (the image dimensions must be less than
128 pixels). Other binary data needs to be converted into a
suitable text representation as well. If you're not sure what to
use, the developers can help you about that point. They will then
also wrap your test case into a Python module for integration
with the framework.

Contributors sufficiently familiar with the Python programming
language and the PyUnit test framework can also submit complete
test suites in Python. Please use the existing tests in the
"py_tests" directory as a template, and check out the helper
modules in "py_tests/unit_tools".

In any case, please note that we can't use any shell scripts or
similar tools in our tests. All tests should be able to run on
all supported platforms, where your favourite shell may not be
available. The Python programming language is available for
pretty much any platform, so we decided to use only that.



