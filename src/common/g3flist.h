/* RCSid $Id: g3flist.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
#ifndef __G3FLIST_H
#define __G3FLIST_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#include "g3vector.h"


#define	G3FL_STARTCAP	50
#define G3FL_ID(fl,id)	(fl->comp_size*id)
#define G3FL_GETID(fl,e) ((e-fl->list)/fl->comp_size)

typedef struct g3FList {
	int			size;
	int			cap;
	int			comp_size;
	int 		ipos;
	g3Float*	list;
} g3FList;

#define 	g3fl_get_size(fl)		(fl->size)
#define		g3fl_get_comp_size(fl)	(fl->comp_size)

int			g3fl_init(g3FList* fl,int cs);
g3FList*	g3fl_create(int cs);
void		g3fl_free(g3FList* fl);

int			g3fl_copy(g3FList* dst,g3FList* src);

void		g3fl_clear(g3FList* fl);

int			g3fl_read(g3FList* fl,FILE* fp);
g3Float*	g3fl_append_new(g3FList* fl);
int			g3fl_append(g3FList* fl,g3Float* el);
g3Float*	g3fl_get(g3FList* fl,int id);
g3Float*	g3fl_get_copy(g3FList* fl,int id);
int			g3fl_remove(g3FList* fl,int id);
int			g3fl_remove_last(g3FList* fl);

g3Float*	g3fl_begin(g3FList* fl);
g3Float*	g3fl_next(g3FList* fl);
void		g3fl_rewind(g3FList* fl);

int			g3fl_sort(g3FList* fl,int cnum);




#ifdef __cplusplus
}
#endif

#endif
