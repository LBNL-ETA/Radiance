#ifndef lint
static const char	RCSid[] = "$Id: sort.c,v 1.5 2003/06/16 14:54:54 greg Exp $";
#endif
/*
 *   Sorting routines for meta-files
 */

#include  "paths.h"
#include  "meta.h"

#ifdef _WIN32
  #include <process.h> /* getpid() */
#endif


#define  PBSIZE  1000		/* max size of sort block */
				/* maxalloc must be >= this */

#define  NFILES  4		/* number of sort files */


static void treemerge(int height, int nt, int nf, PLIST *pl,
		int (*pcmp)(), FILE *ofp);
static void sendsort( PLIST  *pl, int  (*pcmp)());
static void order(PRIMITIVE  *pa[], int  n, PLIST  *pl);
static char * tfname( int  lvl, int num);

/*
 *   This sort routine does a combination quicksort and
 *   n-ary mergesort
 */


void
sort(		/* sort primitives according to pcmp */
FILE  *infp,
int  (*pcmp)()		/* compares pointers to pointers to primitives! */
)
{
 PRIMITIVE  *prims[PBSIZE];		/* pointers to primitives */
 PLIST  primlist;			/* our primitives list */
 int  nprims;
 short  done;

 do  {

    for (nprims = 0; nprims < PBSIZE; nprims++)  {	/* read to global */

       if ((prims[nprims] = palloc()) == NULL)
          error(SYSTEM, "memory exhausted in sort");

       readp(prims[nprims], infp);

       if (isglob(prims[nprims]->com))
	  break;
       }

	qsort(prims, nprims, sizeof(*prims), pcmp);	/* sort pointer array */

	if (nprims < PBSIZE)			/* tack on global if one */
	    nprims++;

	order(prims, nprims, &primlist);	/* make array into list */

	sendsort(&primlist, pcmp);		/* send to merge sorter */

	done = primlist.pbot->com == PEOF;

	plfree(&primlist);			/* free up array */

    }  while (!done);

 }





/*
 *  The following routine merges up to NFILES sorted primitive files
 *     with 0 or 1 primitive lists.  Each set of primitives can end with
 *     0 or 1 global commands.
 */

void
pmergesort(	/* merge sorted files with list */

FILE  *fi[],		/* array of input files */
int  nf,		/* number of input files */
PLIST  *pl,		/* sorted list */
int  (*pcmp)(),		/* comparison function, takes primitive handles */
FILE  *ofp		/* output file */
)
{
    PRIMITIVE  *plp;		/* position in list */
    PRIMITIVE  *pp[NFILES];	/* input primitives */
    int  minf;
    PRIMITIVE  *minp;
    register int i;

    if (pl == NULL)
	plp = NULL;			/* initialize list */
    else
	plp = pl->ptop;

    for (i = 0; i < nf; i++) {		/* initialize input files */
	if ((pp[i] = palloc()) == NULL)
	    error(SYSTEM, "memory exhausted in pmergesort");
	readp(pp[i], fi[i]);
    }

    for ( ; ; ) {

	if (plp != NULL && isprim(plp->com))
	    minp = plp;
	else
	    minp = NULL;

	for (i = 0; i < nf; i++)
	    if (isprim(pp[i]->com) &&
	   		(minp == NULL || (*pcmp)(&pp[i], &minp) < 0))
		minp = pp[minf=i];

	if (minp == NULL)
	    break;

	writep(minp, ofp);

	if (minp == plp)
	    plp = plp->pnext;
	else {
	    fargs(pp[minf]);
	    readp(pp[minf], fi[minf]);
	}
    }

    if (plp != NULL && plp->com != PEOF)
	writep(plp, ofp);

    for (i = 0; i < nf; i++) {
	if (pp[i]->com != PEOF)
	    writep(pp[i], ofp);
	pfree(pp[i]);
    }

}



static void
treemerge(	/* merge into one file */

int  height, int nt, int nf,
PLIST  *pl,
int  (*pcmp)(),
FILE  *ofp
)
{
    FILE  *fi[NFILES], *fp;
    int  i;
    
    if (nf <= NFILES) {

	for (i = 0; i < nf; i++)
	    fi[i] = efopen(tfname(height, i + nt*NFILES), "r");
	
	if ((fp = ofp) == NULL)
	    fp = efopen(tfname(height + 1, nt), "w");
	
	pmergesort(fi, nf, pl, pcmp, fp);
	
	for (i = 0; i < nf; i++) {
	    fclose(fi[i]);
	    unlink(tfname(height, i + nt*NFILES));
	}
	if (ofp == NULL) {
	    writeof(fp);
	    fclose(fp);
	}
	
    } else {
    
        for (i = 0; i < (nf-1)/NFILES; i++)
            treemerge(height, i, NFILES, NULL, pcmp, NULL);
            
	treemerge(height, (nf-1)/NFILES, (nf-1)%NFILES + 1, pl, pcmp, NULL);

	treemerge(height + 1, 0, (nf-1)/NFILES + 1, NULL, pcmp, ofp);
            
    }

}




static void
sendsort(		/* send a sorted list */

PLIST  *pl,
int  (*pcmp)()
)
{
    static int  nf = 0,
    		intree = FALSE;
    FILE  *fp;

    if (isglob(pl->pbot->com)) {

					/* send sorted output to stdout */
	if (intree && nf <= NFILES)

	    treemerge(0, 0, nf, pl, pcmp, stdout);

	else if (nf%NFILES == 0)

	    if (intree) {
		treemerge(0, nf/NFILES - 1, NFILES, pl, pcmp, NULL);
		treemerge(1, 0, nf/NFILES, NULL, pcmp, stdout);
	    } else
		treemerge(1, 0, nf/NFILES, pl, pcmp, stdout);

	else {

	    treemerge(0, nf/NFILES, nf%NFILES, pl, pcmp, NULL);
	    treemerge(1, 0, nf/NFILES + 1, NULL, pcmp, stdout);

	}

	nf = 0;				/* all done */
	intree = FALSE;			/* reset for next time */
	
    } else if (intree && nf%NFILES == 0) {

					/* merge NFILES with list */
	treemerge(0, nf/NFILES - 1, NFILES, pl, pcmp, NULL);
	intree = FALSE;

    } else {

					/* output straight to temp file */
	treemerge(-1, nf++, 0, pl, pcmp, NULL);
	intree = TRUE;

    }



}



static void
order(	/* order the first n array primitives into list */

PRIMITIVE  *pa[],
int  n,
PLIST  *pl
)
{
 register int  i;

 for (i = 1; i < n; i++)
    pa[i-1]->pnext = pa[i];

 pa[n-1]->pnext = NULL;
 pl->ptop = pa[0];
 pl->pbot = pa[n-1];

 }




static char *
tfname(		/* create temporary file name */

int  lvl, int num
)
{
	static char  pathbuf[PATH_MAX];
    static char  fnbuf[PATH_MAX];
	static size_t psiz;

	if (pathbuf[0] == '\0') { /* first time */
		temp_directory(pathbuf, sizeof(pathbuf));
		psiz = strlen(pathbuf);
	}
	snprintf(fnbuf, sizeof(pathbuf)-psiz,
			"%s/S%d%c%d", pathbuf, getpid(), lvl+'A', num);
    /*sprintf(fnbuf, "%sS%d%c%d", TDIR, getpid(), lvl+'A', num);*/

    return(fnbuf);
}
