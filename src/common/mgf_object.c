#ifndef lint
static const char	RCSid[] = "$Id: mgf_object.c,v 3.1 2011/02/18 00:40:25 greg Exp $";
#endif
/*
 * Header file for tracking hierarchical object names
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mgf_parser.h"


int	obj_nnames;			/* depth of name hierarchy */
char	**obj_name;			/* name list */

static int	obj_maxname;		/* allocated list size */

#define ALLOC_INC	16		/* list increment ( > 1 ) */


int
obj_handler(int ac, char **av)		/* handle an object entity statement */
{
	if (ac == 1) {				/* just pop top object */
		if (obj_nnames < 1)
			return(MG_ECNTXT);
		free(obj_name[--obj_nnames]);
		obj_name[obj_nnames] = NULL;
		return(MG_OK);
	}
	if (ac != 2)
		return(MG_EARGC);
	if (!isname(av[1]))
		return(MG_EILL);
	if (obj_nnames >= obj_maxname-1) {	/* enlarge array */
		if (!obj_maxname)
			obj_name = (char **)malloc(
				(obj_maxname=ALLOC_INC)*sizeof(char *));
		else
			obj_name = (char **)realloc(obj_name,
				(obj_maxname+=ALLOC_INC)*sizeof(char *));
		if (obj_name == NULL)
			return(MG_EMEM);
	}
						/* allocate new entry */
	obj_name[obj_nnames] = (char *)malloc(strlen(av[1])+1);
	if (obj_name[obj_nnames] == NULL)
		return(MG_EMEM);
	strcpy(obj_name[obj_nnames++], av[1]);
	obj_name[obj_nnames] = NULL;
	return(MG_OK);
}


void
obj_clear(void)			/* clear object stack */
{
	while (obj_nnames)
		free(obj_name[--obj_nnames]);
	if (obj_maxname) {
		free(obj_name);
		obj_maxname = 0;
	}
}
