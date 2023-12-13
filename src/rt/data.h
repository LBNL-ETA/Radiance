/* RCSid $Id: data.h,v 2.8 2023/12/13 23:26:16 greg Exp $ */
/*
 * Header for data file loading and computation routines.
 */
#ifndef _RAD_DATA_H_
#define _RAD_DATA_H_
#ifdef __cplusplus
extern "C" {
#endif

#define  MAXDDIM	5		/* maximum data dimensions */

#define  DATATYPE	float		/* single precision to save space */
#define  DATATY		'f'		/* format for DATATYPE */
#define  SPECTY		'c'		/* format for SCOLR */

typedef struct datarray {
	char  *name;			/* name of our data */
	short  type;			/* DATATY, SPECTY, RED, GRN or BLU */
	short  nd;			/* number of dimensions */
	struct {
		DATATYPE  org, siz;		/* coordinate domain */
		int  ne;			/* number of elements */
		DATATYPE  *p;			/* point locations */
	} dim[MAXDDIM];			/* dimension specifications */
	union {
		DATATYPE  *d;			/* float data */
		COLR  *c;			/* RGBE data */
		uby8  *s;			/* spectral data */
		void  *p;			/* generic pointer */
	}  arr;				/* the data */
	struct datarray  *next;		/* next array in list */
} DATARRAY;			/* a data array */


extern DATARRAY	*getdata(char *dname);
extern DATARRAY	*getpict(char *pname);
extern DATARRAY *getspec(char *sname);
extern void	freedata(DATARRAY *dta);
extern double	datavalue(DATARRAY *dp, double *pt);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_DATA_H_ */

