#ifndef lint
static const char	RCSid[] = "$Id: holofile.c,v 3.50 2003/07/07 17:21:50 greg Exp $";
#endif
/*
 * Routines for managing holodeck files
 *
 *	9/30/97	GWLarson
 */

#include "copyright.h"

#include <string.h>

#include "holo.h"

#ifndef CACHESIZE
#ifdef SMLMEM
#define CACHESIZE	5
#else
#define CACHESIZE	17	/* default cache size (Mbytes, 0==inf) */
#endif
#endif
#ifndef FREEBEAMS
#define FREEBEAMS	1500	/* maximum beams to free at a time */
#endif
#ifndef PCTFREE
#define PCTFREE		15	/* maximum fraction to free (%) */
#endif
#ifndef MAXFRAGB
#define MAXFRAGB	16	/* fragment blocks/file to track (0==inf) */
#endif
#ifndef FF_DEFAULT
				/* when to free a beam fragment */
#define FF_DEFAULT	(FF_WRITE|FF_KILL)
#endif
#ifndef MINDIRSEL
				/* minimum directory seek length */
#define MINDIRSEL	(4*BUFSIZ/sizeof(BEAMI))
#endif

#ifndef BSD
#define write	writebuf	/* safe i/o routines */
#define read	readbuf
#endif

#define FRAGBLK		512	/* number of fragments to allocate at a time */

int	hdfragflags = FF_DEFAULT;	/* tells when to free fragments */
unsigned	hdcachesize = CACHESIZE*1024*1024;	/* target cache size */
unsigned long	hdclock;		/* clock value */

HOLO	*hdlist[HDMAX+1];	/* holodeck pointers (NULL term.) */

static struct fraglist {
	short	nlinks;		/* number of holodeck sections using us */
	short	writable;	/* 0 read-only, <0 write error encountered */
	int	nfrags;		/* number of known fragments */
	BEAMI	*fi;		/* fragments, descending file position */
	off_t	flen;		/* last known file length */
} *hdfragl;		/* fragment lists, indexed by file descriptor */

static int	nhdfragls;	/* size of hdfragl array */


HOLO *
hdalloc(hproto)		/* allocate and set holodeck section based on grid */
HDGRID	*hproto;
{
	HOLO	hdhead;
	register HOLO	*hp;
	int	n;
				/* copy grid to temporary header */
	memcpy((void *)&hdhead, (void *)hproto, sizeof(HDGRID));
				/* compute grid vectors and sizes */
	hdcompgrid(&hdhead);
				/* allocate header with directory */
	n = sizeof(HOLO)+nbeams(&hdhead)*sizeof(BEAMI);
	if ((hp = (HOLO *)malloc(n)) == NULL)
		return(NULL);
				/* copy header information */
	copystruct(hp, &hdhead);
				/* allocate and clear beam list */
	hp->bl = (BEAM **)malloc((nbeams(hp)+1)*sizeof(BEAM *)+sizeof(BEAM));
	if (hp->bl == NULL) {
		free((void *)hp);
		return(NULL);
	}
	memset((void *)hp->bl, '\0', (nbeams(hp)+1)*sizeof(BEAM *)+sizeof(BEAM));
	hp->bl[0] = (BEAM *)(hp->bl+nbeams(hp)+1);	/* set blglob(hp) */
	hp->fd = -1;
	hp->dirty = 0;
	hp->priv = NULL;
				/* clear beam directory */
	memset((void *)hp->bi, '\0', (nbeams(hp)+1)*sizeof(BEAMI));
	return(hp);		/* all is well */
}


char *
hdrealloc(ptr, siz, rout)	/* (re)allocate memory, retry then error */
char	*ptr;
unsigned	siz;
char	*rout;
{
	register char	*newp;
					/* call malloc/realloc */
	if (ptr == NULL) newp = (char *)malloc(siz);
	else newp = (char *)realloc((void *)ptr, siz);
					/* check success */
	if (newp == NULL && rout != NULL) {
		hdfreecache(25, NULL);	/* free some memory */
		errno = 0;		/* retry */
		newp = hdrealloc(ptr, siz, NULL);
		if (newp == NULL) {	/* give up and report error */
			sprintf(errmsg, "out of memory in %s", rout);
			error(SYSTEM, errmsg);
		}
	}
	return(newp);
}


hdattach(fd, wr)	/* start tracking file fragments for some section */
register int	fd;
int	wr;
{
	if (fd >= nhdfragls) {
		hdfragl = (struct fraglist *)hdrealloc((char *)hdfragl,
				(fd+1)*sizeof(struct fraglist), "hdattach");
		memset((void *)(hdfragl+nhdfragls),
				'\0', (fd+1-nhdfragls)*sizeof(struct fraglist));
		nhdfragls = fd+1;
	}
	hdfragl[fd].nlinks++;
	hdfragl[fd].writable = wr;		/* set writable flag */
	hdfragl[fd].flen = lseek(fd, (off_t)0, 2);	/* get file length */
}


/* Do we need a routine to locate file fragments given known occupants? */


hdrelease(fd)		/* stop tracking file fragments for some section */
register int	fd;
{
	if (fd < 0 | fd >= nhdfragls || !hdfragl[fd].nlinks)
		return;
	if (!--hdfragl[fd].nlinks && hdfragl[fd].nfrags) {
		free((void *)hdfragl[fd].fi);
		hdfragl[fd].fi = NULL;
		hdfragl[fd].nfrags = 0;
	}
}


HOLO *
hdinit(fd, hproto)	/* initialize a holodeck section in a file */
int	fd;			/* corresponding file descriptor */
HDGRID	*hproto;		/* holodeck section grid */
{
	off_t	rtrunc;
	off_t	fpos;
	int	writable;
	register HOLO	*hp;
	register int	n;
					/* prepare for system errors */
	errno = 0;
	if ((fpos = lseek(fd, (off_t)0, 1)) < 0)
		error(SYSTEM, "cannot determine holodeck file position");
	if (hproto == NULL) {		/* assume we're loading it */
		HDGRID	hpr;
						/* load header */
		if (read(fd, (char *)&hpr, sizeof(HDGRID)) != sizeof(HDGRID))
			error(SYSTEM, "cannot load holodeck header");
						/* allocate grid */
		if ((hp = hdalloc(&hpr)) == NULL)
			goto memerr;
						/* load beam directory */
		n = nbeams(hp)*sizeof(BEAMI);
		if (read(fd, (char *)(hp->bi+1), n) != n)
			error(SYSTEM, "failure loading holodeck directory");
						/* check that it's clean */
		for (n = nbeams(hp); n > 0; n--)
			if (hp->bi[n].fo < 0) {
				hp->bi[n].fo = 0;
				error(WARNING, "dirty holodeck section");
				break;
			}
						/* check writability */
		if (fd < nhdfragls && hdfragl[fd].nlinks)
			writable = hdfragl[fd].writable;
		else
			writable = lseek(fd, fpos, 0) == fpos &&
				write(fd, (char *)hp, sizeof(HDGRID)) ==
							sizeof(HDGRID);
	} else {			/* else assume we're creating it */
		if ((hp = hdalloc(hproto)) == NULL)
			goto memerr;
						/* write header and skeleton */
		n = nbeams(hp)*sizeof(BEAMI);
		if (write(fd, (char *)hproto, sizeof(HDGRID)) !=
					sizeof(HDGRID) ||
				write(fd, (char *)(hp->bi+1), n) != n)
			error(SYSTEM, "cannot write header to holodeck file");
		writable = 1;
	}
	hp->fd = fd;	
	hp->dirty = 0;
	biglob(hp)->fo = fpos + sizeof(HDGRID);
					/* start tracking fragments */
	hdattach(fd, writable);
					/* check rays on disk */
	fpos = hdfilen(fd);
	biglob(hp)->nrd = rtrunc = 0;
	for (n = hproto == NULL ? nbeams(hp) : 0; n > 0; n--)
		if (hp->bi[n].nrd)
			if (hp->bi[n].fo+hp->bi[n].nrd*sizeof(RAYVAL) > fpos) {
				rtrunc += hp->bi[n].nrd;
				hp->bi[n].nrd = 0;
			} else
				biglob(hp)->nrd += hp->bi[n].nrd;
	if (rtrunc) {
		sprintf(errmsg, "truncated section, %ld rays lost (%.1f%%)",
				rtrunc, 100.*rtrunc/(rtrunc+biglob(hp)->nrd));
		error(WARNING, errmsg);
	}
					/* add to holodeck list */
	for (n = 0; n < HDMAX; n++)
		if (hdlist[n] == NULL) {
			hdlist[n] = hp;
			break;
		}
					/* all done */
	return(hp);
memerr:
	error(SYSTEM, "cannot allocate holodeck grid");
}


hdmarkdirty(hp, i)		/* mark holodeck directory position dirty */
register HOLO	*hp;
int	i;
{
	static BEAMI	smudge = {0, -1};
	int	mindist, minpos;
	register int	j;

	if (!hp->dirty++) {			/* write smudge first time */
		if (lseek(hp->fd, biglob(hp)->fo+(i-1)*sizeof(BEAMI), 0) < 0
				|| write(hp->fd, (char *)&smudge,
					sizeof(BEAMI)) != sizeof(BEAMI))
			error(SYSTEM, "seek/write error in hdmarkdirty");
		hp->dirseg[0].s = i;
		hp->dirseg[0].n = 1;
		return;
	}
						/* insert into segment list */
	for (j = hp->dirty; j--; ) {
		if (!j || hp->dirseg[j-1].s < i) {
			hp->dirseg[j].s = i;
			hp->dirseg[j].n = 1;
			break;
		}
		copystruct(hp->dirseg+j, hp->dirseg+(j-1));
	}
	do {				/* check neighbors */
		mindist = nbeams(hp);		/* find closest */
		for (j = hp->dirty; --j; )
			if (hp->dirseg[j].s -
					(hp->dirseg[j-1].s + hp->dirseg[j-1].n)
					< mindist) {
				minpos = j;
				mindist = hp->dirseg[j].s -
					(hp->dirseg[j-1].s + hp->dirseg[j-1].n);
			}
						/* good enough? */
		if (hp->dirty <= MAXDIRSE && mindist > MINDIRSEL)
			break;
		j = minpos - 1;			/* coalesce neighbors */
		if (hp->dirseg[j].s + hp->dirseg[j].n <
				hp->dirseg[minpos].s + hp->dirseg[minpos].n)
			hp->dirseg[j].n = hp->dirseg[minpos].s +
					hp->dirseg[minpos].n - hp->dirseg[j].s;
		hp->dirty--;
		while (++j < hp->dirty)		/* close the gap */
			copystruct(hp->dirseg+j, hp->dirseg+(j+1));
	} while (mindist <= MINDIRSEL);
}


int
hdsync(hp, all)			/* update beams and directory on disk */
register HOLO	*hp;
int	all;
{
	register int	j, n;

	if (hp == NULL) {		/* do all holodecks */
		n = 0;
		for (j = 0; hdlist[j] != NULL; j++)
			n += hdsync(hdlist[j], all);
		return(n);
	}
					/* sync the beams */
	for (j = (all ? nbeams(hp) : 0); j > 0; j--)
		if (hp->bl[j] != NULL)
			hdsyncbeam(hp, j);
	if (!hp->dirty)			/* directory clean? */
		return(0);
	errno = 0;			/* write dirty segments */
	for (j = 0; j < hp->dirty; j++) {
		if (lseek(hp->fd, biglob(hp)->fo +
				(hp->dirseg[j].s-1)*sizeof(BEAMI), 0) < 0)
			error(SYSTEM, "cannot seek on holodeck file");
		n = hp->dirseg[j].n * sizeof(BEAMI);
		if (write(hp->fd, (char *)(hp->bi+hp->dirseg[j].s), n) != n)
			error(SYSTEM, "cannot update section directory");
	}
	hp->dirty = 0;			/* all clean */
	return(1);
}


unsigned
hdmemuse(all)		/* return memory usage (in bytes) */
int	all;			/* include overhead (painful) */
{
	long	total = 0;
	register int	i, j;

	for (j = 0; hdlist[j] != NULL; j++) {
		total += blglob(hdlist[j])->nrm * sizeof(RAYVAL);
		if (all) {
			total += sizeof(HOLO) + sizeof(BEAM *) +
					nbeams(hdlist[j]) *
						(sizeof(BEAM *)+sizeof(BEAMI));
			for (i = nbeams(hdlist[j]); i > 0; i--)
				if (hdlist[j]->bl[i] != NULL)
					total += sizeof(BEAM);
		}
	}
	if (all)
		for (j = 0; j < nhdfragls; j++) {
			total += sizeof(struct fraglist);
			if (hdfragl[j].nfrags)
				total += FRAGBLK*sizeof(BEAMI) *
					((hdfragl[j].nfrags-1)/FRAGBLK + 1) ;
		}
	return(total);
}


off_t
hdfilen(fd)		/* return file length for fd */
int	fd;
{
	off_t	fpos, flen;

	if (fd < 0)
		return(-1);
	if (fd >= nhdfragls || !hdfragl[fd].nlinks) {
		if ((fpos = lseek(fd, (off_t)0, 1)) < 0)
			return(-1);
		flen = lseek(fd, (off_t)0, 2);
		lseek(fd, fpos, 0);
		return(flen);
	}
	return(hdfragl[fd].flen);
}


off_t
hdfiluse(fd, all)	/* compute file usage (in bytes) */
int	fd;			/* open file descriptor to check */
int	all;			/* include overhead and unflushed data */
{
	off_t	total = 0;
	register int	i, j;

	for (j = 0; hdlist[j] != NULL; j++) {
		if (hdlist[j]->fd != fd)
			continue;
		total += biglob(hdlist[j])->nrd * sizeof(RAYVAL);
		if (all) {
			for (i = nbeams(hdlist[j]); i > 0; i--)
				if (hdlist[j]->bl[i] != NULL)
					total += sizeof(RAYVAL) *
						(hdlist[j]->bl[i]->nrm -
						hdlist[j]->bi[i].nrd);
			total += sizeof(HDGRID) +
					nbeams(hdlist[j])*sizeof(BEAMI);
		}
	}
	return(total);		/* does not include fragments */
}


RAYVAL *
hdnewrays(hp, i, nr)	/* allocate space for add'l rays and return pointer */
register HOLO	*hp;
register int	i;
int	nr;			/* number of new rays desired */
{
	RAYVAL	*p;
	int	n;

	if (nr <= 0) return(NULL);
	CHECK(i < 1 | i > nbeams(hp),
			CONSISTENCY, "bad beam index given to hdnewrays");
	if (hp->bl[i] != NULL)
		hp->bl[i]->tick = hdclock;	/* preempt swap */
	if (hdcachesize > 0 && hdmemuse(0) >= hdcachesize)
		hdfreecache(PCTFREE, NULL);	/* free some space */
	if (hp->bl[i] == NULL) {		/* allocate (and load) */
		n = hp->bi[i].nrd + nr;
		hp->bl[i] = (BEAM *)hdrealloc(NULL, hdbsiz(n), "hdnewrays");
		blglob(hp)->nrm += n;
		if ((n = hp->bl[i]->nrm = hp->bi[i].nrd)) {
			errno = 0;
			if (lseek(hp->fd, hp->bi[i].fo, 0) < 0)
				error(SYSTEM, "seek error on holodeck file");
			n *= sizeof(RAYVAL);
			if (read(hp->fd, (char *)hdbray(hp->bl[i]), n) != n)
				error(SYSTEM,
				"error reading beam from holodeck file");
		}
	} else {				/* just grow in memory */
		hp->bl[i] = (BEAM *)hdrealloc((char *)hp->bl[i],
				hdbsiz(hp->bl[i]->nrm + nr), "hdnewrays");
		blglob(hp)->nrm += nr;
	}
	if (hdfragflags&FF_ALLOC && hp->bi[i].nrd)
		hdfreefrag(hp, i);		/* relinquish old fragment */
	p = hdbray(hp->bl[i]) + hp->bl[i]->nrm;
	hp->bl[i]->nrm += nr;			/* update in-core structure */
	memset((void *)p, '\0', nr*sizeof(RAYVAL));
	blglob(hp)->tick = hp->bl[i]->tick = hdclock++;	/* update LRU clock */
	return(p);				/* point to new rays */
}


BEAM *
hdgetbeam(hp, i)	/* get beam (from file if necessary) */
register HOLO	*hp;
register int	i;
{
	register int	n;

	CHECK(i < 1 | i > nbeams(hp),
			CONSISTENCY, "bad beam index given to hdgetbeam");
	if (hp->bl[i] == NULL) {		/* load from disk */
		if (!(n = hp->bi[i].nrd))
			return(NULL);
		if (hdcachesize > 0 && hdmemuse(0) >= hdcachesize)
			hdfreecache(PCTFREE, NULL);	/* get free space */
		hp->bl[i] = (BEAM *)hdrealloc(NULL, hdbsiz(n), "hdgetbeam");
		blglob(hp)->nrm += hp->bl[i]->nrm = n;
		errno = 0;
		if (lseek(hp->fd, hp->bi[i].fo, 0) < 0)
			error(SYSTEM, "seek error on holodeck file");
		n *= sizeof(RAYVAL);
		if (read(hp->fd, (char *)hdbray(hp->bl[i]), n) != n)
			error(SYSTEM, "error reading beam from holodeck file");
		if (hdfragflags&FF_READ)
			hdfreefrag(hp, i);	/* relinquish old frag. */
	}
	blglob(hp)->tick = hp->bl[i]->tick = hdclock++;	/* update LRU clock */
	return(hp->bl[i]);
}


int
hdfilord(hb1, hb2)	/* order beams for quick loading */
register HDBEAMI	*hb1, *hb2;
{
	register off_t	c;
				/* residents go first */
	if (hb2->h->bl[hb2->b] != NULL)
		return(hb1->h->bl[hb1->b] == NULL);
	if (hb1->h->bl[hb1->b] != NULL)
		return(-1);
				/* otherwise sort by file descriptor */
	if (hb1->h->fd != hb2->h->fd)
		return(hb1->h->fd - hb2->h->fd);
				/* then by position in file */
	c = hb1->h->bi[hb1->b].fo - hb2->h->bi[hb2->b].fo;
	return(c > 0 ? 1 : c < 0 ? -1 : 0);
}


hdloadbeams(hb, n, bf)	/* load a list of beams in optimal order */
register HDBEAMI	*hb;	/* list gets sorted by hdfilord() */
int	n;			/* list length */
void	(*bf)();		/* callback function (optional) */
{
	unsigned	origcachesize, memuse;
	int	bytesloaded, needbytes, bytes2free;
	register BEAM	*bp;
	register int	i;
					/* precheck consistency */
	if (n <= 0) return;
	for (i = n; i--; )
		if (hb[i].h==NULL || hb[i].b<1 | hb[i].b>nbeams(hb[i].h))
			error(CONSISTENCY, "bad beam in hdloadbeams");
					/* sort list for optimal access */
	qsort((void *)hb, n, sizeof(HDBEAMI), hdfilord);
	bytesloaded = 0;		/* run through loaded beams */
	for ( ; n && (bp = hb->h->bl[hb->b]) != NULL; n--, hb++) {
		bp->tick = hdclock;	/* preempt swap */
		bytesloaded += bp->nrm;
		if (bf != NULL)
			(*bf)(bp, hb);
	}
	bytesloaded *= sizeof(RAYVAL);
	if ((origcachesize = hdcachesize) > 0) {
		needbytes = 0;		/* figure out memory needs */
		for (i = n; i--; )
			needbytes += hb[i].h->bi[hb[i].b].nrd;
		needbytes *= sizeof(RAYVAL);
		do {				/* free enough memory */
			memuse = hdmemuse(0);
			bytes2free = needbytes - (int)(hdcachesize-memuse);
			if (bytes2free > (int)(memuse - bytesloaded))
				bytes2free = memuse - bytesloaded;
		} while (bytes2free > 0 &&
				hdfreecache(100*bytes2free/memuse, NULL) < 0);
		hdcachesize = 0;		/* load beams w/o swap */
	}
	for (i = 0; i < n; i++)
		if ((bp = hdgetbeam(hb[i].h, hb[i].b)) != NULL && bf != NULL)
			(*bf)(bp, hb+i);
	hdcachesize = origcachesize;	/* resume dynamic swapping */
}


int
hdfreefrag(hp, i)			/* free a file fragment */
HOLO	*hp;
int	i;
{
	register BEAMI	*bi = &hp->bi[i];
	register struct fraglist	*f;
	register int	j, k;

	if (bi->nrd <= 0)
		return(0);
	DCHECK(hp->fd < 0 | hp->fd >= nhdfragls || !hdfragl[hp->fd].nlinks,
			CONSISTENCY, "bad file descriptor in hdfreefrag");
	f = &hdfragl[hp->fd];
	if (!f->writable)
		return(0);
	if (f->nfrags % FRAGBLK == 0) {	/* delete empty remnants */
		for (j = k = 0; k < f->nfrags; j++, k++) {
			while (f->fi[k].nrd == 0)
				if (++k >= f->nfrags)
					goto endloop;
			if (k > j)
				copystruct(f->fi+j, f->fi+k);
		}
	endloop:
		f->nfrags = j;
	}
	j = f->nfrags++;		/* allocate a slot in free list */
#if MAXFRAGB
	if (j >= MAXFRAGB*FRAGBLK) {
		f->nfrags = j--;	/* stop list growth */
		if (bi->nrd <= f->fi[j].nrd)
			return(0);	/* new one no better than discard */
	}
#endif
	if (j % FRAGBLK == 0) {		/* more (or less) free list space */
		register BEAMI	*newp;
		if (f->fi == NULL)
			newp = (BEAMI *)malloc((j+FRAGBLK)*sizeof(BEAMI));
		else
			newp = (BEAMI *)realloc((void *)f->fi,
					(j+FRAGBLK)*sizeof(BEAMI));
		if (newp == NULL) {
			f->nfrags--;	/* graceful failure */
			return(0);
		}
		f->fi = newp;
	}
	for ( ; ; j--) {		/* insert in descending list */
		if (!j || bi->fo < f->fi[j-1].fo) {
			f->fi[j].fo = bi->fo;
			f->fi[j].nrd = bi->nrd;
			break;
		}
		copystruct(f->fi+j, f->fi+(j-1));
	}
					/* coalesce adjacent fragments */
						/* successors never empty */
	if (j && f->fi[j-1].fo == f->fi[j].fo + f->fi[j].nrd*sizeof(RAYVAL)) {
		f->fi[j].nrd += f->fi[j-1].nrd;
		f->fi[j-1].nrd = 0;
	}
	for (k = j+1; k < f->nfrags; k++)	/* get non-empty predecessor */
		if (f->fi[k].nrd) {
			if (f->fi[j].fo == f->fi[k].fo +
					f->fi[k].nrd*sizeof(RAYVAL)) {
				f->fi[k].nrd += f->fi[j].nrd;
				f->fi[j].nrd = 0;
			}
			break;
		}
	biglob(hp)->nrd -= bi->nrd;		/* tell fragment it's free */
	bi->nrd = 0;
	bi->fo = 0;
	hdmarkdirty(hp, i);			/* assume we'll reallocate */
	return(1);
}


int
hdfragOK(fd, listlen, listsiz)	/* get fragment list status for file */
int	fd;
int	*listlen;
register int32	*listsiz;
{
	register struct fraglist	*f;
	register int	i;

	if (fd < 0 | fd >= nhdfragls || !(f = &hdfragl[fd])->nlinks)
		return(0);		/* listless */
	if (listlen != NULL)
		*listlen = f->nfrags;
	if (listsiz != NULL)
		for (i = f->nfrags, *listsiz = 0; i--; )
			*listsiz += f->fi[i].nrd;
#if MAXFRAGB
	return(f->nfrags < MAXFRAGB*FRAGBLK);
#else
	return(1);			/* list never fills up */
#endif
}


off_t
hdallocfrag(fd, nrays)		/* allocate a file fragment */
int	fd;
uint32	nrays;
{
	register struct fraglist	*f;
	register int	j;
	off_t	nfo;

	if (nrays == 0)
		return(-1L);
	DCHECK(fd < 0 | fd >= nhdfragls || !hdfragl[fd].nlinks,
			CONSISTENCY, "bad file descriptor in hdallocfrag");
	f = &hdfragl[fd];
	for (j = f->nfrags; j-- > 0; )	/* first fit algorithm */
		if (f->fi[j].nrd >= nrays)
			break;
	if (j < 0) {			/* no fragment -- extend file */
		nfo = f->flen;
		f->flen += nrays*sizeof(RAYVAL);
	} else {			/* else use fragment */
		nfo = f->fi[j].fo;
		f->fi[j].fo += nrays*sizeof(RAYVAL);
		f->fi[j].nrd -= nrays;
	}
	return(nfo);
}


int
hdsyncbeam(hp, i)		/* sync beam in memory with beam on disk */
register HOLO	*hp;
register int	i;
{
	int	fragfreed;
	uint32	nrays;
	unsigned int	n;
	off_t	nfo;
					/* check file status */
	if (hdfragl[hp->fd].writable <= 0)
		return(hdfragl[hp->fd].writable);
	DCHECK(i < 1 | i > nbeams(hp),
			CONSISTENCY, "bad beam index in hdsyncbeam");
					/* is current fragment OK? */
	if (hp->bl[i] == NULL || (nrays = hp->bl[i]->nrm) == hp->bi[i].nrd)
		return(0);
					/* relinquish old fragment? */
	fragfreed = hdfragflags&FF_WRITE && hp->bi[i].nrd && hdfreefrag(hp,i);
	if (nrays) {			/* get and write new fragment */
		nfo = hdallocfrag(hp->fd, nrays);
		errno = 0;
		if (lseek(hp->fd, nfo, 0) < 0)
			error(SYSTEM, "cannot seek on holodeck file");
		n = hp->bl[i]->nrm * sizeof(RAYVAL);
		if (write(hp->fd, (char *)hdbray(hp->bl[i]), n) != n) {
			hdfragl[hp->fd].writable = -1;
			hdsync(NULL, 0);	/* sync directories */
			error(SYSTEM, "write error in hdsyncbeam");
		}
		hp->bi[i].fo = nfo;
	} else
		hp->bi[i].fo = 0;
	biglob(hp)->nrd += nrays - hp->bi[i].nrd;
	hp->bi[i].nrd = nrays;
	if (!fragfreed)
		hdmarkdirty(hp, i);		/* need to flag dir. ent. */
	return(1);
}


int
hdfreebeam(hp, i)		/* free beam, writing if dirty */
register HOLO	*hp;
register int	i;
{
	int	nchanged;

	if (hp == NULL) {		/* clear all holodecks */
		nchanged = 0;
		for (i = 0; hdlist[i] != NULL; i++)
			nchanged += hdfreebeam(hdlist[i], 0);
		return(nchanged);
	}
	if (hdfragl[hp->fd].writable < 0)	/* check for file error */
		return(0);
	if (i == 0) {			/* clear entire holodeck */
		if (blglob(hp)->nrm == 0)
			return(0);		/* already clear */
		nchanged = 0;
		for (i = nbeams(hp); i > 0; i--)
			if (hp->bl[i] != NULL)
				nchanged += hdfreebeam(hp, i);
		DCHECK(blglob(hp)->nrm != 0,
				CONSISTENCY, "bad beam count in hdfreebeam");
		return(nchanged);
	}
	DCHECK(i < 1 | i > nbeams(hp),
			CONSISTENCY, "bad beam index to hdfreebeam");
	if (hp->bl[i] == NULL)
		return(0);
					/* check for additions */
	nchanged = hp->bl[i]->nrm - hp->bi[i].nrd;
	if (nchanged)
		hdsyncbeam(hp, i);		/* write new fragment */
	blglob(hp)->nrm -= hp->bl[i]->nrm;
	free((void *)hp->bl[i]);		/* free memory */
	hp->bl[i] = NULL;
	return(nchanged);
}


int
hdkillbeam(hp, i)		/* delete beam from holodeck */
register HOLO	*hp;
register int	i;
{
	int	nchanged;

	if (hp == NULL) {		/* clobber all holodecks */
		nchanged = 0;
		for (i = 0; hdlist[i] != NULL; i++)
			nchanged += hdkillbeam(hdlist[i], 0);
		return(nchanged);
	}
	if (i == 0) {			/* clobber entire holodeck */
		if (biglob(hp)->nrd == 0 & blglob(hp)->nrm == 0)
			return(0);		/* already empty */
		nchanged = 0;
		nchanged = 0;
		for (i = nbeams(hp); i > 0; i--)
			if (hp->bi[i].nrd > 0 || hp->bl[i] != NULL)
				nchanged += hdkillbeam(hp, i);
		DCHECK(biglob(hp)->nrd != 0 | blglob(hp)->nrm != 0,
				CONSISTENCY, "bad beam count in hdkillbeam");
		return(nchanged);
	}
	DCHECK(i < 1 | i > nbeams(hp), CONSISTENCY,
			"bad beam index to hdkillbeam");
	DCHECK(!hdfragl[hp->fd].writable, CONSISTENCY,
			"hdkillbeam called on read-only holodeck");
	if (hp->bl[i] != NULL) {	/* free memory */
		blglob(hp)->nrm -= nchanged = hp->bl[i]->nrm;
		free((void *)hp->bl[i]);
		hp->bl[i] = NULL;
	} else
		nchanged = hp->bi[i].nrd;
	if (hp->bi[i].nrd && !(hdfragflags&FF_KILL && hdfreefrag(hp,i))) {
		biglob(hp)->nrd -= hp->bi[i].nrd;	/* free failed */
		hp->bi[i].nrd = 0;
		hp->bi[i].fo = 0;
		hdmarkdirty(hp, i);
	}
	return(nchanged);
}


int
hdlrulist(hb, nents, n, hp)	/* add beams from holodeck to LRU list */
register HDBEAMI	*hb;		/* beam list */
int	nents;				/* current list length */
int	n;				/* maximum list length */
register HOLO	*hp;			/* section we're adding from */
{
	register int	i, j;
					/* insert each beam from hp */
	for (i = 1; i <= nbeams(hp); i++) {
		if (hp->bl[i] == NULL)		/* check if loaded */
			continue;
#if 0
		if (hp->bl[i]->tick == hdclock)	/* preempt swap? */
			continue;
#endif
		if ((j = ++nents) >= n)		/* grow list if we can */
			nents--;
		for ( ; ; ) {			/* bubble into place */
			if (!--j || hp->bl[i]->tick >=
					hb[j-1].h->bl[hb[j-1].b]->tick) {
				hb[j].h = hp;
				hb[j].b = i;
				break;
			}
			copystruct(hb+j, hb+(j-1));
		}
	}
	return(nents);			/* return new list length */
}


int
hdfreecache(pct, honly)		/* free up cache space, writing changes */
int	pct;				/* maximum percentage to free */
register HOLO	*honly;			/* NULL means check all */
{
	HDBEAMI	hb[FREEBEAMS];
	int	freetarget;
	int	n;
	register int	i;
#ifdef DEBUG
	unsigned	membefore;

	membefore = hdmemuse(0);
#endif
						/* compute free target */
	freetarget = (honly != NULL) ? blglob(honly)->nrm :
			hdmemuse(0)/sizeof(RAYVAL) ;
	freetarget = freetarget*pct/100;
	if (freetarget <= 0)
		return(0);
						/* find least recently used */
	n = 0;
	if (honly != NULL)
		n = hdlrulist(hb, n, FREEBEAMS, honly);
	else
		for (i = 0; hdlist[i] != NULL; i++)
			n = hdlrulist(hb, n, FREEBEAMS, hdlist[i]);
						/* free LRU beams */
	for (i = 0; i < n; i++) {
		hdfreebeam(hb[i].h, hb[i].b);
		if ((freetarget -= hb[i].h->bi[hb[i].b].nrd) <= 0)
			break;
	}
	hdsync(honly, 0);	/* synchronize directories as necessary */
#ifdef DEBUG
	sprintf(errmsg,
	"%dK before, %dK after hdfreecache (%dK total), %d rays short\n",
		membefore>>10, hdmemuse(0)>>10, hdmemuse(1)>>10, freetarget);
	wputs(errmsg);
#endif
	return(-freetarget);	/* return how far past goal we went */
}


hddone(hp)		/* clean up holodeck section and free */
register HOLO	*hp;		/* NULL means clean up all */
{
	register int	i;

	if (hp == NULL) {		/* NULL means clean up everything */
		while (hdlist[0] != NULL)
			hddone(hdlist[0]);
		free((void *)hdfragl);
		hdfragl = NULL; nhdfragls = 0;
		return;
	}
					/* flush all data and free memory */
	hdflush(hp);
					/* release fragment resources */
	hdrelease(hp->fd);
					/* remove hp from active list */
	for (i = 0; hdlist[i] != NULL; i++)
		if (hdlist[i] == hp) {
			while ((hdlist[i] = hdlist[i+1]) != NULL)
				i++;
			break;
		}
	free((void *)hp->bl);		/* free beam list */
	free((void *)hp);		/* free holodeck struct */
}
