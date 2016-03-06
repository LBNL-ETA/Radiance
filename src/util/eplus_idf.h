/* RCSid $Id: eplus_idf.h,v 2.8 2016/03/06 01:13:18 schorsch Exp $ */
/*
 *  eplus_idf.h
 *
 *  EnergyPlus Input Data File i/o declarations
 *
 *  Created by Greg Ward on 1/31/14.
 */

#ifndef _RAD_EPLUS_IDF_H_
#define _RAD_EPLUS_IDF_H_

#include "platform.h"
#include "lookup.h"

#ifdef __cplusplus
extern "C" {
#endif

#define	IDF_MAXLINE	512			/* maximum line length */
#define IDF_MAXARGL	104			/* maximum argument length */

/* Input Data File object argument list */
typedef struct s_idf_field {
	struct s_idf_field	*next;		/* next in object list */
	char			*rem;		/* string following argument */
	char			val[2];		/* value (extends struct) */
} IDF_FIELD;

/* Input Data File object */
typedef struct s_idf_object {
	const char		*pname;		/* object name (type) */
	struct s_idf_object	*pnext;		/* next object same type */
	struct s_idf_object	*dnext;		/* next object in IDF */
	int			nfield;		/* field count */
	IDF_FIELD		*flist;		/* field list */
	char			rem[1];		/* comment (extends struct) */
} IDF_OBJECT;

/* Input Data File w/ dictionary */
typedef struct {
	char			*hrem;		/* header remarks */
	IDF_OBJECT		*pfirst;	/* first object in file */
	IDF_OBJECT		*plast;		/* last object loaded */
	LUTAB			ptab;		/* object hash table */
} IDF_LOADED;

/* Create a new object with empty field list (comment & prev optional) */
extern IDF_OBJECT	*idf_newobject(IDF_LOADED *idf, const char *pname,
					const char *comm, IDF_OBJECT *prev);

/* Add a field to the given object and follow with the given text */
extern int		idf_addfield(IDF_OBJECT *param, const char *fval,
							const char *comm);

/* Retrieve the indexed field from object (first field index is 1) */
extern IDF_FIELD	*idf_getfield(IDF_OBJECT *param, int fn);

/* Delete the specified object from the IDF */
extern int		idf_delobject(IDF_LOADED *idf, IDF_OBJECT *param);

/* Move the specified object to the given position in the IDF */
extern int		idf_movobject(IDF_LOADED *idf, IDF_OBJECT *param,
						IDF_OBJECT *prev);

/* Get a named object list */
extern IDF_OBJECT	*idf_getobject(IDF_LOADED *idf, const char *pname);

/* Read a object and fields from an open file and add to end of list */
extern IDF_OBJECT	*idf_readobject(IDF_LOADED *idf, FILE *fp);

/* Initialize an IDF struct */
extern IDF_LOADED	*idf_create(const char *hdrcomm);

/* Add comment(s) to header */
extern int		idf_add2hdr(IDF_LOADED *idf, const char *hdrcomm);

/* Load an Input Data File */
extern IDF_LOADED	*idf_load(const char *fname);

/* Write a object and fields to an open file */
int			idf_writeparam(IDF_OBJECT *param, FILE *fp,
					int incl_comm);

/* Write out an Input Data File */
extern int		idf_write(IDF_LOADED *idf, const char *fname,
					int incl_comm);

/* Free a loaded IDF */
extern void		idf_free(IDF_LOADED *idf);

#ifdef __cplusplus
}
#endif
#endif	/* ! _RAD_EPLUS_IDF_H_ */
