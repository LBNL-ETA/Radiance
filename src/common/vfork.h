/* RCSid $Id: vfork.h,v 2.4 2003/02/25 02:47:22 greg Exp $ */
/*
 * Header for routines using vfork() system call.
 */

#include "copyright.h"


#if !defined(BSD) || defined(sparc)

#define	vfork	fork

#endif

extern int	vfork();
