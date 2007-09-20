#ifndef lint
static const char	RCSid[] = "$Id: gencat.c,v 2.6 2007/09/20 21:28:35 greg Exp $";
#endif
/*****************************************************************************
  This program is to make series of right triangles forming hyperbolic cosin
  (ie, cosh) curve in between of 2 points.
       ^             
   f(h)| 
       |		     
       |               0\ pc:(hc, fc)
       |               |    \ 
       |               |       \
       |	       |          \  
       |   pa:(ha, fa) 0-------------0 pb:(hb, fb)
       |
       0------------------------------------------------> h
   
  Given arguments: 
		    material   
		    name
		    (x0, y0, z0), (x1, y1, z1)
		    k        const. value K
		    d        distant length desired between 2 points

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

char  *cmtype, *cname;
double z0, z1;
double k, D;
double d;
double z, h;

#ifdef notdef
double Newton( b)
double b;
{
   if (fabs(cosh(k*D+b)-cosh(b)-(z1-z0)/k) < .001)
      return (b);
   else {
      b = b - (cosh(k*D+b)-cosh(b)-(z1-z0)/k)/(sinh(k*D+b)-sinh(b));
      Newton (b);
   }
}
#endif

double Newton(bl)
double bl;
{
	double b;
	int n = 10000;

	while (n--) {
        	b = bl- (cosh(k*D+bl)-cosh(bl)-(z1-z0)/k)/(sinh(k*D+bl)-sinh(bl));
		if (fabs(b-bl) < .00001)
			return(b);
		bl = b;
	}
fprintf(stderr, "Interation limit exceeded -- invalid K value\n");
	exit(1);
}


static void
printhead(ac, av)		/* print command header */
register int  ac;
register char  **av;
{
	putchar('#');
	while (ac--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');
}

int
main (argc, argv)
int argc;
char *argv[];
{
   double x0, y0;
   double x1, y1;
   double b;
   double delh;
   double f, fprim;
   double hb, hc, fb, fc;
   int  n;

   if (argc != 11) {
      fprintf(stderr, "Usage: %s  material name x0 y0 z0 x1 y1 z1 k d\n",
      		argv[0]);
      exit(1);
   }      

   cmtype = argv[1];
   cname = argv[2];
   x0 = atof(argv[3]); y0 = atof(argv[4]); z0 = atof(argv[5]);
   x1 = atof(argv[6]); y1 = atof(argv[7]); z1 = atof(argv[8]);
   k = atof(argv[9]); d = atof(argv[10]);
   D = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
   b = Newton(0.0); 
   z = z0 - k * cosh(b);
   printhead(argc, argv);

   n = 0;
   for (h=0; h<=D; ) {
      f =  k * cosh(k*h+b) + z;
      fprim = k* k * sinh(k*h+b);
      delh = d / sqrt(1+fprim*fprim); 
      fprim = k * k * sinh(k*(h+delh/2)+b); 
      hb = sqrt(.01*fprim*fprim/(1+fprim*fprim));
      fb = sqrt(.01/(1+(fprim*fprim)));
      hc = sqrt(.04/(1+fprim*fprim));
      fc = sqrt(.04*fprim*fprim/(1+fprim*fprim));

      printf("\n%s polygon %s.%d\n", cmtype, cname, ++n);
      printf("0\n0\n9\n");
      printf("%f %f %f\n", h*(x1-x0)/D+x0, h*(y1-y0)/D+y0, f);
      if (fprim < 0)
      {
         printf("%f %f %f\n", (h+hc)*(x1-x0)/D+x0, (h+hc)*(y1-y0)/D+y0, f-fc);
         printf("%f %f %f\n", (h+hb)*(x1-x0)/D+x0, (h+hb)*(y1-y0)/D+y0, f+fb);
      }
      else if (fprim > 0)
      {
         printf("%f %f %f\n", (h+hc)*(x1-x0)/D+x0, (h+hc)*(y1-y0)/D+y0, f+fc);
         printf("%f %f %f\n", (h-hb)*(x1-x0)/D+x0, (h-hb)*(y1-y0)/D+y0, f+fb);
      }
      else 
      { 
         printf("%f %f %f\n", (h+.2)*(x1-x0)/D+x0, (h+.2)*(y1-y0)/D+y0, f);
         printf("%f %f %f\n", h*(x1-x0)/D+x0, h*(y1-y0)/D+y0, f+.1);
      }
      h += delh;
   }
   return 0;
}

