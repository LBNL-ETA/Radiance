/* RCSid: $Id: rholo.h,v 3.21 2003/02/22 02:07:25 greg Exp $ */
/*
 * Header file for rholo program
 */

#include "holo.h"
#include "vars.h"

#ifndef RPACKSIZ
#define RPACKSIZ	21		/* good packet size */
#endif

typedef struct packet {
	int2	hd;		/* holodeck section (first) */
	int4	bi;		/* beam index (second) */
	int4	nr;		/* number of rays (third) */
	int4	nc;		/* number calculated (fourth) */
	RAYVAL	ra[RPACKSIZ];	/* ray values (fifth) */
	float	*offset;	/* offset array if !vbool(OBSTRUCTIONS) */
	struct packet	*next;	/* next in packet list */
} PACKET;		/* a beam packet */

typedef struct {
	int2	hd;		/* holodeck section (first) */
	int4	bi;		/* beam index (second) */
	int4	nr;		/* number of rays (third) */
	int4	nc;		/* number calculated (fourth) */
} PACKHEAD;		/* followed by ray values */

#define packsiz(nr)	(sizeof(PACKHEAD)+(nr)*sizeof(RAYVAL))
#define packra(p)	((RAYVAL *)((p)+1))

typedef struct {
	int	nb;		/* number of beams in list */
	PACKHEAD	*bl;	/* allocated beam list */
} BEAMLIST;		/* a list of beam requests */

typedef struct {
	FVECT	vpt;		/* view (eye point) position */
	double	rng;		/* desired mean radius for sample rays */
} VIEWPOINT;		/* target eye position */

				/* input variables */
#define CACHE		0		/* amount of memory to use as cache */
#define DISKSPACE	1		/* how much disk space to use */
#define EYESEP		2		/* eye separation distance */
#define GEOMETRY	3		/* section geometry */
#define GRID		4		/* target grid size */
#define OBSTRUCTIONS	5		/* shall we track obstructions? */
#define OCTREE		6		/* octree file name */
#define PORTS		7		/* section portals */
#define RENDER		8		/* rendering options */
#define REPORT		9		/* report interval and error file */
#define RIF		10		/* rad input file */
#define SECTION		11		/* holodeck section boundaries */
#define TIME		12		/* maximum rendering time */
#define VDIST		13		/* virtual distance calculation */

#define NRHVARS		14		/* number of variables */

#define RHVINIT { \
	{"CACHE",	2,	0,	NULL,	fltvalue}, \
	{"DISKSPACE",	3,	0,	NULL,	fltvalue}, \
	{"EYESEP",	3,	0,	NULL,	fltvalue}, \
	{"geometry",	3,	0,	NULL,	NULL}, \
	{"GRID",	2,	0,	NULL,	fltvalue}, \
	{"OBSTRUCTIONS",3,	0,	NULL,	boolvalue}, \
	{"OCTREE",	3,	0,	NULL,	onevalue}, \
	{"portals",	3,	0,	NULL,	NULL}, \
	{"render",	3,	0,	NULL,	catvalues}, \
	{"REPORT",	3,	0,	NULL,	onevalue}, \
	{"RIF",		3,	0,	NULL,	onevalue}, \
	{"section",	3,	0,	NULL,	NULL}, \
	{"TIME",	2,	0,	NULL,	fltvalue}, \
	{"VDISTANCE",	2,	0,	NULL,	boolvalue}, \
}

				/* bundle set requests */
#define BS_NEW		1		/* replace current set with new one */
#define BS_ADD		2		/* add to current set */
#define BS_ADJ		3		/* adjust current set quantities */
#define BS_DEL		4		/* delete from current set */
#define BS_MAX		5		/* set to max of old and new */

extern char	*progname;	/* our program name */
extern char	*hdkfile;	/* holodeck file name */
extern char	froot[];	/* root file name */

extern char	*outdev;	/* output device name */

extern int	readinp;	/* read input from stdin */

extern int	force;		/* allow overwrite of holodeck */

extern int	nowarn;		/* turn warnings off? */

extern int	ncprocs;	/* number of requested compute processes */
extern int	nprocs;		/* number of running compute processes */

extern int	chunkycmp;	/* using "chunky" comparison mode */

extern VIEWPOINT	myeye;	/* target view position */

extern time_t	starttime;	/* time we got started */
extern time_t	endtime;	/* time we should end by */
extern time_t	reporttime;	/* time for next report */

extern long	nraysdone;	/* number of rays done */
extern long	npacksdone;	/* number of packets done */

extern int	rtargc;		/* rtrace command */
extern char	*rtargv[];

extern PACKET	*do_packets(), *get_packets(), *flush_queue();

extern int2	*viewbeams();
