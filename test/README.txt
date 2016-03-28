
Radiance Testing Framework
--------------------------

A toolkit to test all (eventually) components of the Radiance
synthetic image generation system for conformance to their
specification.


Limitations

We use the Python unittest module to run our tests. This means
that we're currently restricted to test only complete programs,
and not actual units (since unittest was designed to test Python
units, not C). A C-level testing framework may be added later.

There's no good way to automatically test GUI programs like
rview. We have to rely on good human testers to check whether
those work correctly or not.


Requirements

You need a working installation of Python 2.7 or 3.x on your
system. Radiance must be either built with the executables still
in the source tree (preferrable to test before installing), or as
a live installation.


How to run tests

The simplest way to run tests is to use the SCons build system.
The file ray/INSTALL.scons explains the requirements and details.
Once you have SCons working, go to the ray directory and type

 $> scons build
 $> scons test

The first command will build Radiance, and place the executables
in a platform-specific directory below ray/scbuild/. 
The second command will automatically execute all available tests
in the environment created by the build.

Other build systems may chose to integrate the tests in a similar
way. The file "run_tests.py" can either be invoked as a script or
imported as a module. Note that in either case, you probably need
to supply the correct paths to the Radiance binaries and library.

As a script:

  usage: run_tests.py [-V] [-H] [-p bindir] [-l radlib] [-c cat]

  optional arguments:
    -V         Verbose: Print all executed test cases to stderr
    -H         Help: print this text to stderr and exit
    -p bindir  Path to Radiance binaries
    -l radlib  Path to Radiance library
    -c cat     Category of tests to run (else all)

As a module:

Call the class run_tests.RadianceTests(...) with suitable arguments:
  bindir=[directory ...] - will be prepended to PATH during tests
  radlib=[directory ...] - will be prepended to RAYPATH during tests
  cat=[category ...]     - only test those categories (else TESTCATS)
  V=False                - if True, verbose listing of executed tests
  
Both methods will run all the tests, or just the category given
as the value of the "cat" argument.


What gets tested

There are several test categories, each containing a number of test
suites, each containing one or more tests. When running tests, each
test category will be printed to the console. Depending on the
settings, the individual test cases may also be listed, or just
indicated with a dot. And last the total results for each category are
shown.

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

The selection of tests to run is still very much incomplete, but
will hopefully grow over time. You can contribute by creating
tests too! Please ask on the code development mailing list first,
so that we can avoid overlaps between the work of different
contributors.

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
use, the developers will be happy to assist you. They will then
also wrap your test case into a Python module for integration
with the framework.

Contributors sufficiently familiar with the Python programming
language and the unittest framework can also submit complete
test suites in Python. Please use the existing tests below the
"testcases" directory as a template, and check out the helper
modules in ".../lib/pyradlib" (where ".../lib" is the location of
the Radiance support library).

Also note the pseudo-builtin module "testsupport" temporarily
created by the RadianceTests() class (see the docstrings there),
which provides information about the various required directory
locations.

And lastly you'll find that we have deliberately included a space
character in the name of the "test data" directory, because it is
a design requirement that all our executables can handle path
names with spaces.

In any case, remember that we can't use any shell scripts or
similar tools in our tests. All tests should be able to run on
all supported platforms, where your favourite shell may not be
available. The Python programming language is available for
pretty much any platform, so we decided to use only that.



