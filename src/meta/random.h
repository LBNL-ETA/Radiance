/* RCSid: $Id: random.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *  random.h - header file for random(3) function.
 *
 *     10/1/85
 */

extern long  random();

#define  frandom()	(random()/2147483648.0)
