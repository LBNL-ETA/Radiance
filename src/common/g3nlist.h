/* RCSid $Id: g3nlist.h,v 2.2 2015/08/18 15:02:53 greg Exp $ */
#ifndef __G3NAMEDLIST_H
#define __G3NAMEDLIST_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include "gbasic.h"
#include "g3list.h"

#define	G3NL_NAME_LEN	200
#define G3NL_STARTCAP	2



typedef	char	g3nl_names[G3NL_NAME_LEN];

typedef	struct	g3nList {
	create_f		ocreate;
	free_f			ofree;
	copy_f			ocopy;
	int				cap;
	int				size;
	unsigned char**			items;
	g3nl_names*		names;
} g3nList;

	

g3nList*	g3nl_create(create_f cf,free_f ff,copy_f cpf);
//void	g3nl_free(g3nList* nl);

//unsigned char*	g3nl_append_new(g3nList* nl,const char* name);
int		g3nl_append(g3nList* nl,const char* name,unsigned char* el);
unsigned char*	g3nl_get(g3nList* nl,int id);
char*			g3nl_get_name(g3nList* nl,int id);
unsigned char*	g3nl_get_by_name(g3nList* nl,const char* name);
//unsigned char*	g3nl_get_copy(g3nList* nl,int id);
//unsigned char*	g3nl_get_copy_by_name(g3nList* nl,const char* name);

int		g3nl_remove(g3nList* nl,int id);
int		g3nl_remove_by_name(g3nList* nl,const char* name);
unsigned char*	g3nl_remove_last(g3nList* nl);
int		g3nl_get_size(g3nList* nl);



#define	G3N_MKHEADER(pref,type) \
	g3nList*    pref ## _nlist_create(); \
	void 	pref ## _nlist_free(g3nList* nl); \
	void 	pref ## _nlist_free_content(g3nList* nl); \
	type*	pref ## _nlist_append_new(g3nList* nl,const char* name); \
	int  	pref ## _nlist_append(g3nList* nl,const char* name,type* it); \
	type*	pref ## _nlist_get(g3nList* nl,int id); \
	char*	pref ## _nlist_get_name(g3nList* nl,int id); \
	type*	pref ## _nlist_get_by_name(g3nList* nl,const char* name); \
	type*	pref ## _nlist_get_copy(g3nList* nl,int id); \
	type*	pref ## _nlist_get_copy_by_name(g3nList* nl,const char* name); \
	int 	pref ## _nlist_remove(g3nList* nl,int id); \
	int		pref ## _nlist_remove_by_name(g3nList* nl,const char* name); \
	type*	pref ## _nlist_remove_last(g3nList* nl); \
	int		pref ## _nlist_get_size(g3nList* nl); \
	void	pref ## _nlist_clear(g3nList* nl); 
	
	

#define	G3N_MKLIST(pref,type,crf,ff,cpf) \
	g3nList*	pref ## _nlist_create() { \
		return g3nl_create((create_f) crf,(free_f) ff,(copy_f) cpf); \
	}; \
	void	pref ## _nlist_free(g3nList* nl) { \
		free(nl->items); \
		free(nl); \
		nl = NULL; \
	} \
	void	pref ## _nlist_free_content(g3nList* nl) { \
		int i; \
		for(i=0;i<nl->cap;i++)\
			nl->ofree((unsigned char*) nl->items[i]);\
		free(nl->items); \
		free(nl); \
		nl = NULL; \
	} \
	type*	pref ## _nlist_append_new(g3nList* nl,const char* name) { \
		type* res = (type*)  nl->ocreate(); \
		if (!g3nl_append(nl,name,((unsigned char*) res))) \
			return NULL; \
		return res; \
	} \
	int		pref ## _nlist_append(g3nList* nl,const char* name,type* it) { \
		return g3nl_append(nl,name,((unsigned char*) it)); \
	} \
	type*	pref ## _nlist_get(g3nList* nl,int id) { \
		return ((type*) g3nl_get(nl,id)); \
	} \
	char*	pref ## _nlist_get_name(g3nList* nl,int id) { \
		return g3nl_get_name(nl,id); \
	} \
	type*	pref ## _nlist_get_by_name(g3nList* nl,const char* name) { \
		return ((type*) g3nl_get_by_name(nl,name)); \
	} \
	type*	pref ## _nlist_get_copy(g3nList* nl,int id) { \
		unsigned char* it; \
		unsigned char* res; \
		if (!(nl->ocopy)) \
			return NULL; \
		res =  nl->ocreate(); \
		if (!(it = (unsigned char*) pref ## _nlist_get(nl,id))) \
			return NULL; \
		if (!(nl->ocopy(res,it))) { \
			fprintf(stderr,"nlist: copying failed for type %s\n","type"); \
			return NULL; \
		} \
		return ((type*) res); \
	} \
	type*	pref ## _nlist_get_copy_by_name(g3nList* nl,const char* name) { \
		int i; \
		unsigned char* it; \
		unsigned char* res; \
		if (!(nl->ocopy)) \
			return NULL; \
		res =  nl->ocreate(); \
		for(i=0;i<nl->size;i++)\
			if (strncmp(nl->names[i],name,G3NL_NAME_LEN-1)) \
				break; \
		if (i == nl->size) \
			return NULL; \
		if (!(it = (unsigned char*) pref ## _nlist_get_copy(nl,i)))  \
				return NULL; \
		if (!(nl->ocopy(res,it))) { \
			fprintf(stderr,"nlist: copying failed for type %s\n","type"); \
			 return NULL; \
		} \
		return ((type*) res); \
	} \
	int	pref ## _nlist_remove(g3nList* nl,int id) { \
		return g3nl_remove(nl,id); \
	}\
	int	pref ## _nlist_remove_by_name(g3nList* nl,const char* name) {\
		return g3nl_remove_by_name(nl,name); \
	} \
	type* pref ## _nlist_remove_last(g3nList* nl) { \
		return ((type*) g3nl_remove_last(nl)); \
	} \
	int	pref ## _nlist_get_size(g3nList* nl) { \
		return g3nl_get_size(nl);  \
	} \
	void pref ## _nlist_clear(g3nList* nl) { \
		nl->size = 0; \
	}



#ifdef __cplusplus
}
#endif
#endif
