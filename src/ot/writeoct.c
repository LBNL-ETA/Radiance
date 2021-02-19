#ifndef lint
static const char RCSid[] = "$Id: writeoct.c,v 2.8 2021/02/19 17:45:37 greg Exp $";
#endif
/*
 *  writeoct.c - routines for writing octree information to stdout.
 *
 *     7/30/85
 */

#include  "standard.h"

#include  "octree.h"
#include  "object.h"
#include  "oconv.h"

#ifdef putc_unlocked		/* avoid horrendous overhead of flockfile */
#undef putc
#define putc    putc_unlocked
#endif

static void oputstr(char *s);
static void putfullnode(OCTREE fn);
static void oputint(long i, int siz);
static void oputflt(double f);
static void puttree(OCTREE ot);


void
writeoct(		/* write octree structures to stdout */
	int  store,
	CUBE  *scene,
	char  *ofn[]
)
{
	char  sbuf[64];
	int  i;
					/* write format number */
	oputint((long)(OCTMAGIC+sizeof(OBJECT)), 2);

	if (!(store & IO_BOUNDS))
		return;
					/* write boundaries */
	for (i = 0; i < 3; i++) {
		sprintf(sbuf, "%.12g", scene->cuorg[i]);
		oputstr(sbuf);
	}
	sprintf(sbuf, "%.12g", scene->cusize);
	oputstr(sbuf);
					/* write object file names */
	if (store & IO_FILES)
		for (i = 0; ofn[i] != NULL; i++)
			oputstr(ofn[i]);
	oputstr("");
					/* write number of objects */
	oputint((long)nobjects, sizeof(OBJECT));

	if (!(store & IO_TREE))
		return;
					/* write the octree */
	puttree(scene->cutree);

	if (fflush(stdout) == EOF)
		error(SYSTEM, "output error in writeoct");

	if (store & IO_FILES || !(store & IO_SCENE))
		return;
					/* write the scene */
	writescene(0, nobjects, stdout);
}


static void
oputstr(			/* write null-terminated string to stdout */
	char  *s
)
{
	if (putstr(s, stdout) == EOF)
		error(SYSTEM, "write error in oputstr");
}


static void
putfullnode(			/* write out a full node */
	OCTREE  fn
)
{
	OBJECT  oset[MAXSET+1];
	int  i;

	objset(oset, fn);
	for (i = 0; i <= oset[0]; i++)
		oputint((long)oset[i], sizeof(OBJECT));
}


static void
oputint(			/* write a siz-byte integer to stdout */
	long  i,
	int  siz
)
{
	if (putint(i, siz, stdout) == EOF)
		error(SYSTEM, "write error in oputint");
}


static void
oputflt(			/* put out floating point number */
	double  f
)
{
	if (putflt(f, stdout) == EOF)
		error(SYSTEM, "write error in oputflt");
}


static void
puttree(			/* write octree to stdout in pre-order form */
	OCTREE  ot
)
{
	int  i;
	
	if (istree(ot)) {
		putc(OT_TREE, stdout);		/* indicate tree */
		for (i = 0; i < 8; i++)		/* write tree */
			puttree(octkid(ot, i));
	} else if (isfull(ot)) {
		putc(OT_FULL, stdout);		/* indicate fullnode */
		putfullnode(ot);		/* write fullnode */
	} else
		putc(OT_EMPTY, stdout);		/* indicate empty */
}
