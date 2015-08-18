/* RCSid $Id: maxheap.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef	__MAXHEAP_H
#define __MAXHEAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "g3vector.h"

#define	HP_STARTCAP	5000001

#define	HP_INVALID	0
#define	HP_BUILD	1
#define	HP_SORTED	2

typedef	struct	heapElem {
	g3Float		key;
	void*		entity;
} heapElem;

typedef	struct maxHeap {
	int 		size;
	int			cap;
	int			state;
	heapElem*	heap;
} maxHeap;

#define	mheap_get_max(hp) 	(hp->heap[0])
#define	mh_left(i)			(2*i + 1)
#define	mh_right(i)			(2*i + 2)
#define	mh_parent(i)		(((int) ((i + 1)/2)) - 1)

int			mheap_init(maxHeap* hp);
void		mheap_free(maxHeap* hp);
maxHeap*	mheap_create();
void		mheap_build(maxHeap* hp);
void		mheap_sort(maxHeap* hp);

heapElem	mheap_remove_max(maxHeap* hp);
int			mheap_insert(maxHeap* hp,g3Float key,void* entity);


#ifdef __cplusplus
}
#endif

#endif
