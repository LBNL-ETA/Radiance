/*
 *  Disk2Square.c
 *  
 * Code by Peter Shirley and Kenneth Chiu
 *
 * Associated paper in ~/Documents/Others/Shirley+Chiu_JGT-97.pdf
 *
 * Modified interface slightly (G. Ward)
 */

#include <math.h>

/*
 This transforms points on [0,1]^2 to points on unit disk centered at
 origin.  Each "pie-slice" quadrant of square is handled as a seperate
 case.  The bad floating point cases are all handled appropriately.
 The regions for (a,b) are:

        phi = pi/2
       -----*-----
       |\       /|
       |  \ 2 /  |
       |   \ /   |
 phi=pi* 3  *  1 *phi = 0
       |   / \   |
       |  / 4 \  |
       |/       \|
       -----*-----
        phi = 3pi/2

change log:
    26 feb 2004.  bug fix in computation of phi (now this matches the paper,
                  which is correct).  thanks to Atilim Cetin for catching this.
*/


/* Map a [0,1]^2 square to a unit radius disk */
void
SDsquare2disk(double ds[2], double seedx, double seedy)
{

   double phi, r;
   double a = 2.*seedx - 1;   /* (a,b) is now on [-1,1]^2 */
   double b = 2.*seedy - 1;

   if (a > -b) {     /* region 1 or 2 */
       if (a > b) {  /* region 1, also |a| > |b| */
           r = a;
           phi = (M_PI/4.) * (b/a);
       }
       else       {  /* region 2, also |b| > |a| */
           r = b;
           phi = (M_PI/4.) * (2. - (a/b));
       }
   }
   else {        /* region 3 or 4 */
       if (a < b) {  /* region 3, also |a| >= |b|, a != 0 */
            r = -a;
            phi = (M_PI/4.) * (4. + (b/a));
       }
       else       {  /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
            r = -b;
            if (b != 0.)
                phi = (M_PI/4.) * (6. - (a/b));
            else
                phi = 0.;
       }
   }

   ds[0] = r * cos(phi);
   ds[1] = r * sin(phi);

}

/* Map point on unit disk to a unit square in [0,1]^2 range */
void
SDdisk2square(double sq[2], double diskx, double disky)
{
    double r = sqrt( diskx*diskx + disky*disky );
    double phi = atan2( disky, diskx );
    double a, b;

    if (phi < -M_PI/4) phi += 2*M_PI;  /* in range [-pi/4,7pi/4] */
    if ( phi < M_PI/4) {         /* region 1 */
        a = r;
        b = phi * a / (M_PI/4);
    }
    else if ( phi < 3*M_PI/4 ) { /* region 2 */
        b = r;
        a = -(phi - M_PI/2) * b / (M_PI/4);
    }
    else if ( phi < 5*M_PI/4 ) { /* region 3 */
        a = -r;
        b = (phi - M_PI) * a / (M_PI/4);
    }
    else {                       /* region 4 */
        b = -r;
        a = -(phi - 3*M_PI/2) * b / (M_PI/4);
    }

    sq[0] = (a + 1) * 0.5;
    sq[1] = (b + 1) * 0.5;
}
