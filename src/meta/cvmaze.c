#ifndef lint
static const char	RCSid[] = "$Id: cvmaze.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *      THE MAZE GAME  definitions and global structures
 *
 *           Written by Greg Ward 8/4/84
 */


#include  "meta.h"


#define  ERROR  (-1)		/* something wrong */

#define  OK  0			/* everything swell */

#define  TRUE  1

#define  FALSE  0

#define  UP  1

#define  DOWN  2

#define  LEFT  3

#define  RIGHT  4

#define  min(a, b)  ((a) < (b) ? (a) : (b))

#define  max(a, b)  ((a) > (b) ? (a) : (b))

#define  VSEP  '|'		/* vertical separator */

#define  HSEP  '-'		/* horizontal separator */

#define  ISECT  '+'		/* intersection character */

#define  MARKS  " ^v<>*"	/* maze markings */

#define  TUP  1			/* mark for trace up */

#define  TDOWN  2		/* mark for trace down */

#define  TLEFT  3		/* mark for trace left */

#define  TRIGHT  4		/* mark for trace right */

#define  MAN  5			/* mark for player */

#define  MAXSIDE  50		/* maximum maze side */

#define  YES  1			/* edge is there */

#define  NO  0			/* edge isn't there */

#define  MAYBE  2		/* undecided */





#define  DIRECTION  char

#define  EDGE  char

#define  MRKINDX  char

struct  MAZE  {
               EDGE  r0[MAXSIDE+1][MAXSIDE],
                     c0[MAXSIDE+1][MAXSIDE];

               int  nrows, ncols;

               MRKINDX  mark[MAXSIDE][MAXSIDE];

               int  rt, rb, cl, cr;

               };



/*
 *  Convert maze to metafile
 */





struct MAZE  amaze;	/* our working maze */

int  cellsize;		/* the size of a maze cell (in metacoordinates) */

char  *progname;

#define  CPMEOF  0x1a


main(argc, argv)

int  argc;
char  **argv;

{
    FILE  *fp;

#ifdef  CPM
    fixargs("cvmaze", &argc, &argv);
#endif
    progname = *argv++;
    argc--;
        
    if (argc)
    
        for ( ; argc; argc--, argv++) {
            if ((fp = fopen(*argv, "r")) == NULL) {
                fprintf(stderr, "%s:  no such file\n", *argv);
                exit(1);
            }
            execute(fp);
            fclose(fp);
        }
        
    else
    
        execute(stdin);

    writeof(stdout);
    
    return(0);
}




execute(fp)		/* process file fp */

FILE  *fp;

{
    
    while (readmaze(fp, &amaze, MARKS) != ERROR) {
   
       cellsize = (XYSIZE-1)/max(amaze.nrows, amaze.ncols);
       printmaze(stdout, &amaze, MARKS);
       pglob(PEOP, 0200, NULL);
    }
    
}





printmaze(fp, mzp, mrks)		/* print maze to fp */

FILE  *fp;
struct MAZE  *mzp;
char  *mrks;

{
 int  i;
 
 for (i = mzp->rt; i < mzp->rb; i++)  {
 
    putrline(fp, mzp, i);
    
    putcline(fp, mzp, i, mrks);
    
    }

 putrline(fp, mzp, mzp->rb);
 
 }
 
 
 
 
putrline(fp, mzp, row)			/* print horizontal row to fp */

FILE  *fp;
struct MAZE  *mzp;
int  row;

{
 register int  j;

 for (j = mzp->cl; j < mzp->cr; j++) 
   if (mzp->r0[row][j] == YES)
       plseg(0, j*cellsize, (XYSIZE-1)-row*cellsize, (j+1)*cellsize, (XYSIZE-1)-row*cellsize);
   
 }
 
 
 
 
 
putcline(fp, mzp, row, mrks)		/* print column line to fp */

FILE  *fp;
struct MAZE  *mzp;
int  row;
char  mrks[];

{
 register int  j;
 char  s[2];
 
 s[1] = '\0';
 
 for (j = mzp->cl; j <= mzp->cr; j++)
    if (mzp->c0[j][row] == YES)
       plseg(0, j*cellsize, (XYSIZE-1)-row*cellsize, j*cellsize, (XYSIZE-1)-(row+1)*cellsize);

 }
 
 
 
    

readmaze(fp, mzp, mrks)			/* read maze from fp into mzp */

FILE  *fp;
struct MAZE  *mzp;
char  *mrks;

{
 int  i;

 mzp->rt = mzp->rb = mzp->cl = mzp->cr = mzp->nrows = mzp ->ncols = 0;

 if ((mzp->cr = mzp->ncols = getrline(fp, mzp, 0)) == 0)
    return(ERROR);

 for (i = 0; i <= MAXSIDE; i++)  {

    if (getcline(fp, mzp, i, mrks) != mzp->ncols)
       break;

    if (getrline(fp, mzp, i+1) != mzp->ncols)
       return(ERROR);

    }

 if (i > MAXSIDE)
    return(ERROR);

 mzp->rb = mzp->nrows = i;
 return(OK);
 }
 


getrline(fp, mzp, row)			/* get row from fp */

FILE  *fp;
struct MAZE  *mzp;
int  row;

{
 int  j, nc;
 char  *fgets(), linbuf[4 + 4*MAXSIDE];

 if (fgets(linbuf, sizeof(linbuf), fp) == NULL)
    return(0);

 nc = (strlen(linbuf) - 2) / 4;

 for (j = 0; j < nc; j++)

    switch (linbuf[2 + 4*j])  {

       case HSEP:
          mzp->r0[row][j] = YES;
          break;

       case ' ':
          mzp->r0[row][j] = NO;
          break;

       default:
          return(0);
          break;

       }

 return(nc);
 }





getcline(fp, mzp, row, mrks)		/* get column line from fp */

FILE  *fp;
struct MAZE  *mzp;
int  row;
char  mrks[];

{
 int  j, nc, k;
 char  *fgets(), linbuf[4 + 4*MAXSIDE];

 if (fgets(linbuf, sizeof(linbuf), fp) == NULL)
    return(0);

 nc = (strlen(linbuf) - 2) / 4;

 for (j = 0; j < nc; j++)  {

    switch (linbuf[4*j])  {

       case VSEP:
          mzp->c0[j][row] = YES;
          break;

       case ' ':
          mzp->c0[j][row] = NO;
          break;

       default:
          return(0);
          break;

       }

    for (k = 0; mrks[k]; k++)

       if (linbuf[2 + 4*j] == mrks[k])  {
          mzp->mark[row][j] = k;
          break;
          }

    if (!mrks[k])
       return(0);

    }

 switch (linbuf[4*nc])  {

    case VSEP:
       mzp->c0[nc][row] = YES;
       break;

    case ' ':
       mzp->c0[nc][row] = NO;
       break;

    default:
       return(0);
       break;

    }

 return(nc);
 }
