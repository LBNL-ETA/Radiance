#ifndef lint
static const char	RCSid[] = "$Id: sm_usets.c,v 3.2 2003/02/22 02:07:25 greg Exp $";
#endif
/*
 * Quadtree-specific set operations with unsorted sets.
 */

#include "standard.h"
#include "sm_flag.h"
#include "object.h"
#include "sm_qtree.h"


#define QTSETIBLK	2048     	/* set index allocation block size */
#define QTELEMBLK	8		/* set element allocation block */
#define QTELEMMOD(ne)	((ne)&7)	/* (ne) % QTELEMBLK */

#define QTONTHRESH(ne)	(!QTELEMMOD((ne)+1))
#define QTNODESIZ(ne)	((ne) + QTELEMBLK - QTELEMMOD(ne))

OBJECT	**qtsettab= NULL;		/* quadtree leaf node table */
QUADTREE  qtnumsets=0;			/* number of used set indices */
static int  qtfreesets = EMPTY;		/* free set index list */
int4 *qtsetflag = NULL;
static int qtallocsets =0;

qtclearsetflags()
{

  if(!qtsetflag)
    return;

    bzero((char *)qtsetflag,FLAG_BYTES(qtallocsets));
}


QUADTREE
qtnewleaf(oset)			/* return new leaf node for object set */
OBJECT  *oset;
{
	register QUADTREE  osi;

	if (*oset <= 0)
		return(EMPTY);		/* should be error? */
	if (qtfreesets != EMPTY) {
		osi = qtfreesets;
		qtfreesets = (int)qtsettab[osi];
	} else if ((osi = qtnumsets++) % QTSETIBLK == 0) {
		qtsettab = (OBJECT **)realloc((char *)qtsettab,
				(unsigned)(osi+QTSETIBLK)*sizeof(OBJECT *));
		if (qtsettab == NULL)
			goto memerr;
		qtsetflag = (int4 *)realloc((char *)qtsetflag,
				   FLAG_BYTES(osi+ QTSETIBLK));
		if (qtsetflag == NULL)
			goto memerr;
		if(qtallocsets)
		  bzero((char *)((char *)(qtsetflag)+FLAG_BYTES(qtallocsets)),
		      FLAG_BYTES(osi+QTSETIBLK)-FLAG_BYTES(qtallocsets));
		else
		  bzero((char *)(qtsetflag),FLAG_BYTES(osi +QTSETIBLK));

		qtallocsets = osi + QTSETIBLK; 
	}
	qtsettab[osi] = (OBJECT *)malloc(QTNODESIZ(*oset)*sizeof(OBJECT));
	if (qtsettab[osi] == NULL)
		goto memerr;
	setcopy(qtsettab[osi], oset);
	return(QT_INDEX(osi));
memerr:
	error(SYSTEM, "out of memory in qtnewleaf\n");
}


deletuelem(os, obj)		/* delete obj from unsorted os, no questions */
register OBJECT  *os;
OBJECT  obj;
{
	register int  i,j;
	OBJECT *optr;
	/* Take 2nd to last and put in position of obj: move last down:
	   want to preserve last added
	 */
	i = (*os)--;
#ifdef DEBUG
	if( i <= 0)
	  error(CONSISTENCY,"deleteuelem():empty set\n");
#endif
	if(i < 3)
	{
	  if((i == 2) && (*(++os) == obj))
	    *os = *(os+1);
	  return;
	}
	optr = os + i;
	os++;
	while (i-- > 1 && *os != obj)
	   os++;
#ifdef DEBUG
	if( *os != obj)
	  error(CONSISTENCY,"deleteuelem():id not found\n");
#endif
	*os = *(optr-1);
	*(optr-1) = *optr;
}

QUADTREE
qtdeletuelem(qt, id)		/* delete element from unsorted leaf node */
QUADTREE  qt;
OBJECT  id;
{
	register QUADTREE  lf;
#ifdef DEBUG
	if(id < 0)
		eputs("qtdeletuelem():Invalid triangle id\n");
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "empty/bad node in qtdelelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	if (qtsettab[lf][0] <= 1) {		/* blow leaf away */
		free((void *)qtsettab[lf]);
		qtsettab[lf] = (OBJECT *)qtfreesets;
		qtfreesets = lf;
		return(EMPTY);
	}
	deletuelem(qtsettab[lf], id);
	if (QTONTHRESH(qtsettab[lf][0]))
		qtsettab[lf] = (OBJECT *)realloc((char *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0])*sizeof(OBJECT));
	return(qt);
}

QUADTREE
qtdelelem(qt, id)		/* delete element from leaf node */
QUADTREE  qt;
OBJECT  id;
{
	register QUADTREE  lf;
#ifdef DEBUG
	if(id < 0)
		eputs("qtdelelem():Invalid triangle id\n");
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "empty/bad node in qtdelelem");
#else
	lf = QT_SET_INDEX(qt);
#endif

	if (!inset(qtsettab[lf], id))
	{
#ifdef DEBUG
	  eputs("id not in leaf in qtdelelem\n");
#endif
	  return(qt);
	}
	if (qtsettab[lf][0] <= 1) {		/* blow leaf away */
		free((void *)qtsettab[lf]);
		qtsettab[lf] = (OBJECT *)qtfreesets;
		qtfreesets = lf;
		return(EMPTY);
	}
	deletelem(qtsettab[lf], id);
	if (QTONTHRESH(qtsettab[lf][0]))
		qtsettab[lf] = (OBJECT *)realloc((char *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0])*sizeof(OBJECT));
	return(qt);
}


QUADTREE
qtaddelem(qt, id)		/* add element to leaf node */
QUADTREE  qt;
OBJECT  id;
{
    OBJECT	newset[2];
    register QUADTREE  lf;

#ifdef DEBUG
	if(id < 0)
		eputs("qtaddelem():Invalid triangle id\n");
#endif
	if (QT_IS_EMPTY(qt)) {		/* create new leaf */
		newset[0] = 1; newset[1] = id;
		return(qtnewleaf(newset));
	}				/* else add element */
#ifdef DEBUG
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtaddelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	
	if (inset(qtsettab[lf], id))
	{
#ifdef DEBUG
	  eputs("id already in leaf in qtaddelem\n");
#endif
	  return(qt);
	}
	if (QTONTHRESH(qtsettab[lf][0])) {
		qtsettab[lf] = (OBJECT *)realloc((char *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0]+1)*sizeof(OBJECT));
		if (qtsettab[lf] == NULL)
			error(SYSTEM, "out of memory in qtaddelem");
	}
	insertelem(qtsettab[lf], id);
	return(qt);
}

#define insertuelem(os,id) ((os)[++((os)[0])] = (id))

int
qtcompressuelem(qt,compress_set)
  QUADTREE qt;
  int (*compress_set)();
{
    register QUADTREE  lf;
    int i,j,osize;
    OBJECT *os;

    if(QT_IS_EMPTY(qt)) 
      return(qt);
#ifdef DEBUG
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtaddelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	os = qtsettab[lf];
	osize = os[0];
	if((i=compress_set(os)) < osize)
	{
	  qtsettab[lf] = (OBJECT *)realloc((char *)qtsettab[lf],
					   QTNODESIZ(i+1)*sizeof(OBJECT));
	  if (qtsettab[lf] == NULL)
	    error(SYSTEM, "out of memory in qtaddelem");
	}
	return(qt);
}


QUADTREE
qtadduelem(qt, id)		/* add element to unsorted leaf node */
QUADTREE  qt;
OBJECT  id;
{
    OBJECT	newset[2];
    register QUADTREE  lf;

#ifdef DEBUG
	if(id < 0)
		eputs("qtadduelem():Invalid sample id\n");
#endif
	if (QT_IS_EMPTY(qt)) {		/* create new leaf */
		newset[0] = 1; newset[1] = id;
		return(qtnewleaf(newset));
	}				/* else add element */
#ifdef DEBUG
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtaddelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	if (QTONTHRESH(qtsettab[lf][0])) {
		qtsettab[lf] = (OBJECT *)realloc((char *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0]+1)*sizeof(OBJECT));
		if (qtsettab[lf] == NULL)
			error(SYSTEM, "out of memory in qtaddelem");
	}
	insertuelem(qtsettab[lf], id);
	return(qt);
}


#ifdef DEBUG
OBJECT *
qtqueryset(qt)			/* return object set for leaf node */
QUADTREE  qt;
{
	register QUADTREE  lf;

	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtqueryset");
	return(qtsettab[lf]);
}
#endif

int
qtinuset(qt,id)
QUADTREE qt;
OBJECT id;
{
  OBJECT *os;
  int i;

  os = qtqueryset(qt);
  for(i=os[0],os++; i > 0; i--,os++)
    if(*os==id)
      return(1);

  return(0);
}


qtfreeleaf(qt)			/* free set and leaf node */
QUADTREE  qt;
{
	register QUADTREE	osi;

	if (!QT_IS_LEAF(qt))
		return;
	osi = QT_SET_INDEX(qt);
	if (osi >= qtnumsets)
		return;
	free((void *)qtsettab[osi]);
	qtsettab[osi] = (OBJECT *)qtfreesets;
	qtfreesets = osi;
}



qtfreeleaves()			/* free ALL sets and leaf nodes */
{
	register int	i;

	while ((i = qtfreesets) != EMPTY) {
		qtfreesets = (int)qtsettab[i];
		qtsettab[i] = NULL;
	}
	for (i = qtnumsets; i--; )
		if (qtsettab[i] != NULL)
			free((void *)qtsettab[i]);
	free((void *)qtsettab);
	qtsettab = NULL;
	if(qtsetflag)
	{
	  free((void *)qtsetflag);
	  qtsetflag=0;
	}
	qtnumsets = 0;
}


/* no bounds on os csptr. This routine is conservative: the result returned
in os is the set intersection, and therefore is bounded by the size of os.
cs is bound to be <= QT_MAXCSET
*/ 
check_set_large(os, cs)  /* modify checked set and set to check */
     register OBJECT  *os;                   /* os' = os - cs */
     register OBJECT  *cs;                   /* cs' = cs + os */
{
    OBJECT  cset[QT_MAXCSET+1];
    register int  i, j;
    int  k;
                                        /* copy os in place, cset <- cs */

    cset[0] = 0;
    k = 0;
    for (i = j = 1; i <= os[0]; i++) {
	while (j <= cs[0] && cs[j] < os[i])
	  if(cset[0] < QT_MAXCSET)
	    cset[++cset[0]] = cs[j++];
	  else
	    j++;
	if (j > cs[0] || os[i] != cs[j]) {      /* object to check */
	    os[++k] = os[i];
	    if(cset[0] < QT_MAXCSET)
	      cset[++cset[0]] = os[i];
	}
    }
    if (!(os[0] = k))               /* new "to check" set size */
       return;                 /* special case */

    while (j <= cs[0])              /* get the rest of cs */
      if(cset[0] == QT_MAXCSET)
	break;
      else
	cset[++cset[0]] = cs[j++];
    
    os = cset;                 /* copy cset back to cs */
    for (i = os[0]; i-- >= 0; )
       *cs++ = *os++;  

}




