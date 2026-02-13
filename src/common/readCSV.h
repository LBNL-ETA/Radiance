/* RCSid $Id$ */
/* Declarations and data structures for readCSV.c */

#ifndef _READCSV_H_
#define _READCSV_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_csvrec {
	struct s_csvrec		*next;		/* next record in list */
	int			nf;		/* field count */
	char			*f[1];		/* field list (extends struct) */
} CSVREC;

extern CSVREC	*alloc_csvrec(int len);
extern CSVREC	*realloc_csvrec(CSVREC *rp, int newlen);
extern int	set_csvfield(CSVREC *rp, int i, const char *fstr);
extern void	free_csv(CSVREC *rp);
extern CSVREC	*read_csvrec(FILE *fp, CSVREC *toappend);
extern CSVREC	*read_csvfile(char *fname);

#define	get_csvfield(rp,i)	(const char *)( (rp!=NULL) & ((i)>=0) && \
					(i)<(rp)->nf ? (rp)->f[i] : NULL )

#ifdef __cplusplus
}
#endif
#endif	/* ! _READCSV_H_ */
