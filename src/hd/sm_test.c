#ifndef lint
static const char	RCSid[] = "$Id: sm_test.c,v 3.1 2003/02/22 02:07:25 greg Exp $";
#endif
/*
 * sm_test.c
 * 
 * A simple OGL display program to test sm code
 */
#include "standard.h"

#include <glut.h>
#include "ui.h"
#include "sm_geom.h"
#include "sm_qtree.h"
#include "sm_stree.h"
#include "sm.h"
#include "sm_draw.h"

struct driver odev;
char *progname;
double
test_normalize(v)			/* normalize a vector, return old magnitude */
register FVECT  v;
{
	register double  len;
	
	len = DOT(v, v);
	
	if (len <= 0.0)
		return(0.0);
	
	if (len <= 1.0+FTINY && len >= 1.0-FTINY)
		len = 0.5 + 0.5*len;	/* first order approximation */
	else
		len = sqrt(len);

	v[0] /= len;
	v[1] /= len;
	v[2] /= len;

	return(len);
}

void
print_usage(char *program)
{
    fprintf(stderr, "usage: %s \n", program);
}


/* Parse the command line arguments to get the input file name. The depth
 * thresholding size may also be specified, as well as the disparity
 * threshold, specified in number of depth ranges, for creating depth
 * disparity lines.
 */
void
parse_cmd_line(int argc, char* argv[])
{
  int i;
  if(argc < 2)
    return;
  
  for(i=1;i<argc;i++)
  if(argv[i][0] == '-')
    switch(argv[i][1]) {
    case 't':
      Tri_intersect = TRUE;
      break;
    case 'd':
      Pt_delete = TRUE;
      break;
    case 'f':
      if(i+1 < argc)
      {
	Filename = &(argv[i+1][0]);
	i++;
      }
      else
	fprintf(stderr,"Must include file name\n");
    }
    /* Nothing for now */
}

void
display()
{
    glutMainLoop();
}

VIEW View = {0,{0,0,0},{0,0,-1},{0,1,0},60,60,0};
int
main(int argc, char *argv[])
{
    COLR colr;
    FVECT v0,v1,v2,a,b,c,dir;
    double d[3];
    double phi,theta;

    progname = argv[0];
    parse_cmd_line(argc,argv);    
    /* Open a display window */
    smDraw_open_window(argc,argv);
   /* First read in the command line arguments */

    odev.v = View;

   if(Filename)
    {
      next_sample(a,a,colr,1);
      /* Initialize the data structures for 10 samples*/
      smInit(262000);
    }
    else
      smInit(100);
    display();
}













