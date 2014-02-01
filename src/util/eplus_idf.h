/* RCSid $Id: eplus_idf.h,v 2.1 2014/02/01 01:28:43 greg Exp $ */
/*
 *  eplus_idf.h
 *
 *  EnergyPlus Input Data File i/o declarations
 *
 *  Created by Greg Ward on 1/31/14.
 */

#ifndef _RAD_EPLUS_IDF_H_
#define _RAD_EPLUS_IDF_H_

#include "lookup.h"

#ifdef __cplusplus
extern "C" {
#endif

#define	IDF_MAXLINE	512			/* maximum line length */
#define IDF_MAXARGL	128			/* maximum argument length */

/* Input Data File parameter argument list */
typedef struct s_idf_field {
	struct s_idf_field	*next;		/* next in parameter list */
	char			*rem;		/* string following argument */
	char			arg[2];		/* argument (extends struct) */
} IDF_FIELD;

/* Input Data File parameter */
typedef struct s_idf_parameter {
	const char		*pname;		/* parameter name (type) */
	struct s_idf_parameter	*pnext;		/* next parameter same type */
	struct s_idf_parameter	*dnext;		/* next parameter in IDF */
	IDF_FIELD		*flist;		/* field list */
	char			rem[1];		/* comment (extends struct) */
} IDF_PARAMETER;

/* Input Data File w/ dictionary */
typedef struct {
	char			*hrem;		/* header remarks */
	IDF_PARAMETER		*pfirst;	/* first parameter in file */
	IDF_PARAMETER		*plast;		/* last parameter loaded */
	LUTAB			ptab;		/* parameter table */
} IDF_LOADED;

/* Create a new parameter with empty field list (comment optional) */
extern IDF_PARAMETER	*idf_newparam(IDF_LOADED *idf, const char *pname,
							const char *comm,
						IDF_PARAMETER *prev);

/* Add a field to the given parameter and follow with the given text */
extern int		idf_addfield(IDF_PARAMETER *param, const char *fval,
							const char *comm);

/* Delete the specified parameter from our IDF */
extern int		idf_delparam(IDF_LOADED *idf, IDF_PARAMETER *param);

/* Get a named parameter list */
extern IDF_PARAMETER	*idf_getparam(IDF_LOADED *idf, const char *pname);

/* Read a parameter and fields from an open file and add to end of list */
extern IDF_PARAMETER	*idf_readparam(IDF_LOADED *idf, FILE *fp);

/* Initialize an IDF struct */
extern IDF_LOADED	*idf_create(const char *hdrcomm);

/* Load an Input Data File */
extern IDF_LOADED	*idf_load(const char *fname);

/* Write a parameter and fields to an open file */
int			idf_writeparam(IDF_PARAMETER *param, FILE *fp);

/* Write out an Input Data File */
extern int		idf_write(IDF_LOADED *idf, const char *fname);

/* Free a loaded IDF */
extern void		idf_free(IDF_LOADED *idf);

#ifdef __cplusplus
}
#endif
#endif	/* ! _RAD_EPLUS_IDF_H_ */
