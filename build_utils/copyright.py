
import sys
import string

_ltext = '''

                     RADIANCE LICENSE AGREEMENT

Radiance is a registered copyright of The Regents of the University of
California ("The Regents"). The Regents grant to you a nonexclusive,
nontransferable license ("License") to use Radiance source code without
fee.  You may not sell or distribute Radiance to others without the
prior express written permission of The Regents.  You may compile and
use this software on any machines to which you have personal access,
and may share its use with others who have access to the same machines.

NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR ANY
OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS
OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.  By downloading, using
or copying this software, you agree to abide by the intellectual property laws
and all other applicable laws of the United States, and by the terms of this
License Agreement. Ownership of the software shall remain solely in The
Regents.  The Regents shall have the right to terminate this License
immediately by written notice upon your breach of, or noncompliance with, any
of its terms.  You shall be liable for any infringement or damages resulting
from your failure to abide by the terms of this License Agreement.



'''

def show_license():
	sys.stderr.write(_ltext)
	sys.stderr.write(
	'Do you understand and accept the terms of this agreement [n]?\n\n')

	answer = ''
	s =  'Please enter "yes" or "no", or use ^C to exit: '
	while answer not in ['y', 'ye', 'yes', 'n', 'no']:
		if answer: sys.stderr.write('invalid input "%s"\n' % answer)
		answer = string.lower(raw_input(s))

	if answer[0] == 'y':
		return
	else:
		sys.stderr.write('\n*** Installation cancelled ***\n')
		sys.exit(1)

