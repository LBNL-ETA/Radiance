#ifndef lint
static const char RCSid[] = "$Id: maxheap.c,v 2.3 2018/04/19 15:31:27 greg Exp $";
#endif
#include <stdio.h>
#include <stdlib.h>

#include "maxheap.h"


static	int	resize_heap(maxHeap* hp,int sz)
{
	int i;

	if (sz == hp->cap)
		return 1;
	hp->heap = (heapElem*) realloc(hp->heap,sz*sizeof(heapElem));
	if (!hp->heap) 
		return 0;
	for(i=hp->size;i<sz;i++) {
		hp->heap[i].entity = NULL;
		hp->heap[i].key = -FHUGE;
	}
	hp->cap = sz;
	if (sz < hp->size) {
		hp->size = sz;
		hp->state = HP_INVALID;
	}
	if (sz == 0) {
		hp->heap = NULL;
		hp->state = HP_BUILD;
	}
	return 1;
}

int		mheap_init(maxHeap* hp)
{
	hp->size = 0;
	hp->cap	 = 0;
	hp->heap = NULL;
	hp->state = HP_BUILD;
	return resize_heap(hp,HP_STARTCAP);
}

void	mheap_free(maxHeap* hp)
{
	resize_heap(hp,0);
	free(hp);
}

maxHeap*	mheap_create()
{
	maxHeap* res = (maxHeap*) calloc(1,sizeof(maxHeap));
	if (!res) {
		fprintf(stderr,"Error allocating heap struct\n");
		return NULL;
	}
	if (!mheap_init(res))
		return NULL;
	return res;
}
		

static	void	swap_elem(heapElem* e1,heapElem* e2)
{
	heapElem e;

	e = *e1;
	*e1 = *e2;
	*e2 = e;
}

static	void	heapify(maxHeap* hp,int i)
{
	int	l,r,m;

	l = mh_left(i);
	r = mh_right(i);
	if ((l < hp->size) && (hp->heap[l].key > hp->heap[i].key)) 
		m = l;
	else
		m = i;
	if ((r < hp->size) && (hp->heap[m].key < hp->heap[r].key))
		m = r;
	if (m != i) {
		swap_elem(&(hp->heap[m]),&(hp->heap[i]));
		heapify(hp,m);
	}
}

void	mheap_build(maxHeap* hp)
{
	int i;

	for(i=((int) hp->size/2)-1;i>=0;i--) 
		heapify(hp,i);
	hp->state = HP_BUILD;
}

void	mheap_sort(maxHeap* hp)
{
	int i,sz;

	if (hp->state == HP_SORTED) 
		return;
	if (hp->state != HP_BUILD)
		mheap_build(hp);
	sz = hp->size;
	for(i=(hp->size-1);i>0;i--) {
		swap_elem(&(hp->heap[0]),&(hp->heap[i]));
		hp->size--;
		heapify(hp,0);
	}
	hp->size = sz;
	hp->state = HP_SORTED;
}
		
		
heapElem	mheap_remove_max(maxHeap* hp)
{
	heapElem res;
	if (hp->size < 1) {
		res.key = -FHUGE;
		res.entity = NULL;
		return res;
	}	
	res = hp->heap[0];
	hp->heap[0] = hp->heap[hp->size - 1];
	hp->size--;
	heapify(hp,0);
	return res;
}

static	int	inc_key(maxHeap* hp,int i,g3Float k)
{
	if (k < hp->heap[i].key) {
		fprintf(stderr,"New key error in inc_key(mheap)\n");
		hp->state = HP_INVALID;
		return 0;
	}
	hp->heap[i].key = k;
	while((i > 0) && (hp->heap[mh_parent(i)].key < hp->heap[i].key)) {
		swap_elem(&(hp->heap[mh_parent(i)]),&(hp->heap[i]));
		i = mh_parent(i);
	}
	return 1;
}

int		mheap_insert(maxHeap* hp,g3Float key,void* entity)
{
	if (hp->size == hp->cap) {
		if (!resize_heap(hp,hp->cap*2)) {
			hp->state = HP_INVALID;
			return 0;
		}
	}
	hp->heap[hp->size].entity = entity;
	hp->heap[hp->size].key = -FHUGE;
	hp->size++;
	return inc_key(hp,(hp->size - 1),key);
}


#ifdef	TEST_MHEAP
#include "timeutil.h"
		
int		main(int argc,char** argv)
{
	g3Float* k;
	int i,n;
	maxHeap* hp;
	heapElem e;
	double t;
	hp = mheap_create();

	if (argc > 1)
		n = atoi(argv[1]);
	else
		n = 10;
	
	t = get_time();
	for(i=0;i<n;i++) {
		k = (g3Float*) malloc(sizeof(g3Float));
		*k = 10.0*rand()/RAND_MAX;
		if (!mheap_insert(hp,*k,k)) {
			fprintf(stderr,"Error inserting\n");
			return EXIT_FAILURE;
		}
	}
	t = get_time() - t;
	fprintf(stderr,"insert: %d %f\n",n,t);
	t = get_time();
	mheap_sort(hp);
	t = get_time() - t;
	fprintf(stderr,"sort: %d %f\n",n,t);
	for(i=0;i<n;i++) 
		printf("%d %f \n",i,hp->heap[i].key);
	printf("-------------\n");
	t = get_time();
	mheap_build(hp);
	t = get_time() - t;
	fprintf(stderr,"build: %d %f\n",n,t);
	for(i=0;i<n;i++) {
		e = mheap_remove_max(hp);
		printf("%d %f \n",hp->size,e.key);
		if ((i < (n-1)) && (e.key < mheap_get_max(hp).key)) {
			fprintf(stderr,"oooops: not sorted\n");
		}
	}
	return EXIT_SUCCESS;
}

#endif
