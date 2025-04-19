#ifndef lint
static const char RCSid[] = "$Id: cmconst.c,v 2.1 2025/04/19 17:12:59 greg Exp $";
#endif
/*
 * Constants referenced from both cmatrix.o and rmatrix.o
 */

#include "color.h"

const char	stdin_name[] = "<stdin>";

const char	*cm_fmt_id[] = {
			"unknown", COLRFMT, CIEFMT, SPECFMT,
			"float", "ascii", "double"
		};

const int	cm_elem_size[] = {
			0, 4, 4, 0, 3*sizeof(float), 0, 3*sizeof(double)
		};
