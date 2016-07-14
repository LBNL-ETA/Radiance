#ifndef lint
static const char RCSid[] = "$Id: g3flist.c,v 2.3 2016/07/14 17:32:12 greg Exp $";
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "g3flist.h"


static int	resize_list(g3FList* fl,int sz)
{
	fl->list = (g3Float*) realloc(fl->list,(sz*fl->comp_size*sizeof(g3Float)));
	if (!fl->list)
		return 0;
	fl->cap = sz;
	if (sz < fl->size)
		fl->size = sz;
	if (sz == 0)
		fl->list = NULL;
	return 1;
}

int		g3fl_init(g3FList* fl,int cs)
{
	fl->cap = 0;
	fl->size = 0;
	fl->list = NULL;
	fl->comp_size = cs;
	fl->ipos = 0;
	return resize_list(fl,G3FL_STARTCAP);
}

g3FList*	g3fl_create(int cs)
{
	g3FList* res = (g3FList*) malloc(sizeof(g3FList));
	if (!g3fl_init(res,cs))
		return NULL;
	return res;
}

void	g3fl_free(g3FList* fl)
{
	free(fl->list);
	free(fl);
}

g3Float*	g3fl_append_new(g3FList* fl)
{
	if (fl->size == fl->cap) {
		if (!resize_list(fl,fl->cap*2)) {
			fprintf(stderr,"Out of memory in g3fl_append_new\n");
			return NULL;
		}
	}
	fl->size++;
	return (&fl->list[G3FL_ID(fl,(fl->size-1))]);
}

int		g3fl_append(g3FList* fl,g3Float* el)
{
	g3Float* ne;

	if (!(ne = g3fl_append_new(fl)))
		return 0;
	memcpy(ne,el,(fl->comp_size*sizeof(g3Float)));
	return 1;
}

g3Float*	g3fl_get(g3FList* fl,int id)
{
	return &(fl->list[G3FL_ID(fl,id)]);
}

g3Float*	g3fl_begin(g3FList* fl)
{
	if (fl->size == 0)
		return NULL;
	fl->ipos = 1;
	return g3fl_get(fl,0);
}

void	g3fl_rewind(g3FList* fl)
{
	fl->ipos = 0;
}

static int comp_sort(const void* __v1,const void* __v2) {
	const g3Float *v1 = __v1;
	const g3Float* v2 = __v2;
	if (v1[0] < v2[0])
		return -1;
	if (v1[0] > v2[0])
		return 1;
	return 0;
}

int         g3fl_sort(g3FList* fl,int cnum)
{
	g3Float* n;
    if (cnum >= fl->comp_size)
        return 0;
	if (cnum > 0) /* brute force method to avoid using inline functions */
		for(n=g3fl_begin(fl);n != NULL;n=g3fl_next(fl))
			gb_swap(n[0], n[cnum]);
    qsort(fl->list, fl->size, fl->comp_size*sizeof(g3Float), comp_sort);
	if (cnum > 0) /* brute force method to avoid using inline functions */
		for(n=g3fl_begin(fl);n != NULL;n=g3fl_next(fl))
			gb_swap(n[0], n[cnum]);
    return 1;
}


g3Float*	g3fl_next(g3FList* fl)
{
	if (fl->ipos >= fl->size)
		return NULL;
	return g3fl_get(fl,fl->ipos++);
}

g3Float*	g3fl_get_copy(g3FList* fl,int id)
{
	g3Float* ne;

	if (!(ne = (g3Float*) malloc(fl->comp_size*sizeof(g3Float))))
		return NULL;
	memcpy(ne,&(fl->list[G3FL_ID(fl,id)]),(fl->comp_size*sizeof(g3Float)));
	return ne;
}

int		g3fl_remove(g3FList* fl,int id)
{
	if ((id < 0) || (id >= fl->size))
		return 0;
	if (id < (fl->size - 1)) 
		memmove(&(fl->list[G3FL_ID(fl,id)]),&(fl->list[G3FL_ID(fl,id+1)]),
				fl->comp_size*sizeof(g3Float));
	fl->size--;
	return 1;
}

int		g3fl_remove_last(g3FList* fl)
{
	fl->size--;
	if (fl->size < (fl->cap/2))
		return resize_list(fl,(fl->cap/2));
	return 1;
}


#ifdef TEST_FLIST

int main(int argc,char** argv)
{
	int i;
	g3FList* fl = g3fl_create(2);
	g3Float* n;
	
	for(i=0;i<10;i++) {
		n = g3fl_append_new(fl);
		n[0] = i;
		n[1] = 10-i;
	}
	g3fl_sort(fl, 1);
	for(n=g3fl_begin(fl);n != NULL;n=g3fl_next(fl))
		printf("%f %f\n",n[0],n[1]);
	g3fl_free(fl);
	return EXIT_SUCCESS;
} 
		

#endif
