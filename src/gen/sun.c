#ifndef lint
static const char	RCSid[] = "$Id: sun.c,v 2.5 2009/06/06 20:18:32 greg Exp $";
#endif
/*
 *           SOLAR CALCULATIONS
 *
 *               3/31/87
 *
 */

#include  <math.h>

#ifdef M_PI
#define  PI	M_PI
#else
#define  PI	3.141592653589793
#endif

double  s_latitude = 0.66;	/* site latitude (radians) */
double  s_longitude = 2.13;	/* site longitude (radians) */
double  s_meridian = 2.0944;	/* standard meridian (radians) */


int
jdate(		/* Julian date (days into year) */
	int month,
	int day
)
{
	static short  mo_da[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
	
	return(mo_da[month-1] + day);
}


double
stadj(		/* solar time adjustment from Julian date */
	int  jd
)
{
	return( 0.170 * sin( (4*PI/373) * (jd - 80) ) -
		0.129 * sin( (2*PI/355) * (jd - 8) ) +
		12 * (s_meridian - s_longitude) / PI );
}


double
sdec(		/* solar declination angle from Julian date */
	int  jd
)
{
	return( 0.4093 * sin( (2*PI/368) * (jd - 81) ) );
}


double
salt(	/* solar altitude from solar declination and solar time */
	double sd,
	double st
)
{
	return( asin( sin(s_latitude) * sin(sd) -
			cos(s_latitude) * cos(sd) * cos(st*(PI/12)) ) );
}


double
sazi(	/* solar azimuth from solar declination and solar time */
	double sd,
	double st
)
{
	return( -atan2( cos(sd)*sin(st*(PI/12)),
 			-cos(s_latitude)*sin(sd) -
 			sin(s_latitude)*cos(sd)*cos(st*(PI/12)) ) );
}
