/* RCSid $Id: vfork.h,v 2.5 2003/06/27 06:53:22 greg Exp $ */
/*
 * Header for routines using vfork() system call.
 */

#if !defined(BSD) || defined(sparc)

#define	vfork	fork

#endif

extern int	vfork();
