#ifndef lint
static const char	RCSid[] = "$Id: noise3.c,v 2.14 2018/11/21 20:23:15 greg Exp $";
#endif
/*
 *  noise3.c - noise functions for random textures.
 *
 *     Credit for the smooth algorithm goes to Ken Perlin, and code
 *	translation/implementation to Rahul Narain.
 *	(ref. Improving Noise, Computer Graphics; Vol. 35 No. 3., 2002)
 */

#include "copyright.h"

#include  "ray.h"
#include  "func.h"

static char  noise_name[4][8] = {"noise3x", "noise3y", "noise3z", "noise3"};
static char  fnoise_name[] = "fnoise3";

#define  EPSILON	.0005		/* error allowed in fractal */

#define  frand3(x,y,z)	frand(17*(x)+23*(y)+29*(z))

static double l_noise3(char *nam);
static double noise3(double xnew[3], int i);
static double noise3partial(double f3, double x[3], int i);
static double perlin_noise (double x, double y, double z);
static double frand(long s);
static double fnoise3(double x[3]);


static double
l_noise3(			/* compute a noise function */
	char  *nam
)
{
	int  i;
	double  x[3];
					/* get point */
	x[0] = argument(1);
	x[1] = argument(2);
	x[2] = argument(3);
					/* make appropriate call */
	if (nam == fnoise_name)
		return(fnoise3(x));
	i = 4;
	while (i--)
		if (nam == noise_name[i])
			return(noise3(x,i));
	eputs(nam);
	eputs(": called l_noise3!\n");
	quit(1);
	return 1; /* pro forma return */
}


void
setnoisefuncs(void)			/* add noise functions to library */
{
	int  i;

	funset(fnoise_name, 3, ':', l_noise3);
	i = 4;
	while (i--)
		funset(noise_name[i], 3, ':', l_noise3);
}


static double
frand(			/* get random number from seed */
	long  s
)
{
	s = s<<13 ^ s;
	return(1.0-((s*(s*s*15731+789221)+1376312589)&0x7fffffff)/1073741824.0);
}


static double
fnoise3(			/* compute fractal noise function */
	double  x[3]
)
{
	long  t[3], v[3], beg[3];
	double  fval[8], fc;
	int  branch;
	long  s;
	int  i, j;
						/* get starting cube */
	s = (long)(1.0/EPSILON);
	for (i = 0; i < 3; i++) {
		t[i] = s*x[i];
		beg[i] = s*floor(x[i]);
	}
	for (j = 0; j < 8; j++) {
		for (i = 0; i < 3; i++) {
			v[i] = beg[i];
			if (j & 1<<i)
				v[i] += s;
		}
		fval[j] = frand3(v[0],v[1],v[2]);
	}
						/* compute fractal */
	for ( ; ; ) {
		fc = 0.0;
		for (j = 0; j < 8; j++)
			fc += fval[j];
		fc *= 0.125;
		if ((s >>= 1) == 0)
			return(fc);		/* close enough */
		branch = 0;
		for (i = 0; i < 3; i++) {	/* do center */
			v[i] = beg[i] + s;
			if (t[i] > v[i]) {
				branch |= 1<<i;
			}
		}
		fc += s*EPSILON*frand3(v[0],v[1],v[2]);
		fval[~branch & 7] = fc;
		for (i = 0; i < 3; i++) {	/* do faces */
			if (branch & 1<<i)
				v[i] += s;
			else
				v[i] -= s;
			fc = 0.0;
			for (j = 0; j < 8; j++)
				if (~(j^branch) & 1<<i)
					fc += fval[j];
			fc = 0.25*fc + s*EPSILON*frand3(v[0],v[1],v[2]);
			fval[~(branch^1<<i) & 7] = fc;
			v[i] = beg[i] + s;
		}
		for (i = 0; i < 3; i++) {	/* do edges */
			if ((j = i+1) == 3) j = 0;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			if (++j == 3) j = 0;
			if (branch & 1<<j)
				v[j] += s;
			else
				v[j] -= s;
			fc = fval[branch & ~(1<<i)];
			fc += fval[branch | 1<<i];
			fc = 0.5*fc + s*EPSILON*frand3(v[0],v[1],v[2]);
			fval[branch^1<<i] = fc;
			if ((j = i+1) == 3) j = 0;
			v[j] = beg[j] + s;
			if (++j == 3) j = 0;
			v[j] = beg[j] + s;
		}
		for (i = 0; i < 3; i++)		/* new cube */
			if (branch & 1<<i)
				beg[i] += s;
	}
}


static double
noise3(			/* compute the revised Perlin noise function */
	double  xnew[3], int i
)
{
	static int     gotV;
	static double  x[3];
	static double  f[4];

	if (!gotV || xnew[0] != x[0] || (xnew[1] != x[1]) | (xnew[2] != x[2])) {
		f[3] = perlin_noise(x[0]=xnew[0], x[1]=xnew[1], x[2]=xnew[2]);
		gotV = 0x8;
	}
	if (!(gotV>>i & 1)) {
		f[i] = noise3partial(f[3], x, i);
		gotV |= 1<<i;
	}
	return(f[i]);
}

static double
noise3partial(		/* compute partial derivative for ith coordinate */
	double f3, double x[3], int i
)
{
	double	fc;

	switch (i) {
	case 0:
		fc = perlin_noise(x[0]-EPSILON, x[1], x[2]);
		break;
	case 1:
		fc = perlin_noise(x[0], x[1]-EPSILON, x[2]);
		break;
	case 2:
		fc = perlin_noise(x[0], x[1], x[2]-EPSILON);
		break;
	default:
		return(.0);
	}
	return((f3 - fc)/EPSILON);
}

#define fade(t) ((t)*(t)*(t)*((t)*((t)*6. - 15.) + 10.))

static double lerp(double t, double a, double b) {return a + t*(b - a);}

static double
grad(int hash, double x, double y, double z)
{
    int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
    double u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
                 v = h<4 ? y : h==12|h==14 ? x : z;
    return (!(h&1) ? u : -u) + (!(h&2) ? v : -v);
}

static const int permutation[256] = {151,160,137,91,90,15,
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

#define	p(i)	permutation[(i)&0xff]

static double
perlin_noise(double x, double y, double z) 
{
    int		X, Y, Z;
    double	u, v, w;
    int		A, AA, AB, B, BA, BB;

    X = (int)x-(x<0);                    // FIND UNIT CUBE THAT
    Y = (int)y-(y<0);                    // CONTAINS POINT.
    Z = (int)z-(z<0);
    x -= (double)X;                          // FIND RELATIVE X,Y,Z
    y -= (double)Y;                          // OF POINT IN CUBE.
    z -= (double)Z;
    X &= 0xff; Y &= 0xff; Z &= 0xff;

    u = fade(x);                                  // COMPUTE FADE CURVES
    v = fade(y);                                  // FOR EACH OF X,Y,Z.
    w = fade(z);

    A = p(X  )+Y; AA = p(A)+Z; AB = p(A+1)+Z;          // HASH COORDINATES OF
    B = p(X+1)+Y; BA = p(B)+Z; BB = p(B+1)+Z;          // THE 8 CUBE CORNERS,

    return lerp(w, lerp(v, lerp(u, grad(p(AA  ), x  , y  , z  ),  // AND ADD
                                   grad(p(BA  ), x-1, y  , z  )), // BLENDED
                           lerp(u, grad(p(AB  ), x  , y-1, z  ),  // RESULTS
                                   grad(p(BB  ), x-1, y-1, z  ))),// FROM  8
                   lerp(v, lerp(u, grad(p(AA+1), x  , y  , z-1),  // CORNERS
                                   grad(p(BA+1), x-1, y  , z-1)), // OF CUBE
                           lerp(u, grad(p(AB+1), x  , y-1, z-1),
                                   grad(p(BB+1), x-1, y-1, z-1))));
}
