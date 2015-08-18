/* RCSid $Id: g3list.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
#ifndef __G3LIST_H
#define __G3LIST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "gbasic.h"

#define G3L_STARTCAP	2


typedef	unsigned char*	(*create_f)();
typedef	void	(*free_f)(unsigned char*);
typedef	int		(*copy_f)(unsigned char*,unsigned char*);

typedef	struct	g3List {
	create_f		ocreate;
	free_f			ofree;
	copy_f			ocopy;
	int				cap;
	int				size;
	unsigned char**			items;
} g3List;

	

g3List*	g3l_create(create_f cf,free_f ff,copy_f cpf);

int		g3l_append(g3List* nl,unsigned char* el);
unsigned char*	g3l_get(g3List* nl,int id);

int		g3l_remove(g3List* nl,int id);
unsigned char*	g3l_remove_last(g3List* nl);
int		g3l_get_size(g3List* nl);
void	g3l_invert(g3List* nl);



#define	G3L_MKHEADER(pref,type) \
	g3List*    pref ## _list_create(); \
	void 	pref ## _list_free(g3List* nl); \
	void 	pref ## _list_free_content(g3List* nl); \
	type*	pref ## _list_append_new(g3List* nl); \
	int  	pref ## _list_append(g3List* nl,type* it); \
	type*	pref ## _list_get(g3List* nl,int id); \
	type*	pref ## _list_get_copy(g3List* nl,int id); \
	int		pref ## _list_copy(g3List* dst,g3List* src); \
	int 	pref ## _list_remove(g3List* nl,int id); \
	type*	pref ## _list_remove_last(g3List* nl); \
	int		pref ## _list_get_size(g3List* nl); \
	void	pref ## _list_clear(g3List* nl); \
	int		pref ## _list_index(g3List* nl,type*); \
	int		pref ## _list_replace(g3List* nl,int,type*); 
	
	

#define	G3L_MKLIST(pref,type,crf,ff,cpf) \
	g3List*	pref ## _list_create() { \
		return g3l_create((create_f) crf,(free_f) ff,(copy_f) cpf); \
	}; \
	void	pref ## _list_free(g3List* nl) { \
		free(nl->items); \
		free(nl); \
		nl = NULL; \
	} \
	void	pref ## _list_free_content(g3List* nl) { \
		int i; \
		for(i=0;i<nl->cap;i++)\
			nl->ofree((unsigned char*) nl->items[i]);\
		free(nl->items); \
		free(nl); \
		nl = NULL; \
	} \
	type*	pref ## _list_append_new(g3List* nl) { \
		type* res; \
		res =  (type*)  nl->ocreate(); \
		if (!g3l_append(nl,((unsigned char*) res))) \
			return NULL; \
		return res; \
	} \
	int		pref ## _list_append(g3List* nl,type* it) { \
		return g3l_append(nl,((unsigned char*) it)); \
	} \
	type*	pref ## _list_get(g3List* nl,int id) { \
		return ((type*) g3l_get(nl,id)); \
	} \
	type*	pref ## _list_get_copy(g3List* nl,int id) { \
		unsigned char* it; \
		unsigned char* res; \
		if (!(nl->ocopy)) \
			return NULL; \
		res =  nl->ocreate(); \
		if (!(it = (unsigned char*) pref ## _list_get(nl,id))) \
			return NULL; \
		if (!(nl->ocopy(res,it))) { \
			fprintf(stderr,"list: copying failed for type %s\n","type"); \
			return NULL; \
		} \
		return ((type*) res); \
	} \
	int pref ## _list_copy(g3List* dst,g3List* src) {\
		int i; \
		dst->size = 0; \
		for(i=0;i<src->size;i++)\
			pref ## _list_append(dst,pref ## _list_get_copy(src,i)); \
		return 1; \
	}\
	int	pref ## _list_remove(g3List* nl,int id) { \
		return g3l_remove(nl,id); \
	}\
	type* pref ## _list_remove_last(g3List* nl) { \
		return ((type*) g3l_remove_last(nl)); \
	} \
	int	pref ## _list_get_size(g3List* nl) { \
		return g3l_get_size(nl);  \
	} \
	void pref ## _list_clear(g3List* nl) { \
		nl->size = 0; \
	} \
	int		pref ## _list_index(g3List* nl,type* el) { \
		int i; \
		for(i=0;i<nl->size;i++)\
			if (el == pref ## _list_get(nl,i)) \
				return i; \
		return -1; \
	} \
	int		pref ## _list_replace(g3List* nl,int id,type* el) { \
		if ((id < 0) || (id >= nl->size)) \
			return 0; \
		nl->items[id] = (unsigned char*) el; \
		return 1; \
	}




#ifdef __cplusplus
}
#endif
#endif
