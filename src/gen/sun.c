#ifndef lint
static const char	RCSid[] = "$Id: sun.c,v 2.7 2019/11/07 23:15:07 greg Exp $";
#endif
/*
 *           SOLAR CALCULATIONS
 *
 *               3/31/87
 *
 *		Michalsky algorithm added October 2019:
 *		"The Astronomical Almanac's Algorithm for Approximate
 *		Solar Position (1950-2050)" by Joseph J. Michalsky,
 *		published in 1988 in Solar Energy, Vol. 40, No. 3.
 *		Also added correction to sdec() (365 was 368 originally)
 *
 */

#include  <math.h>
#include  "sun.h"

#ifdef M_PI
#define  PI	M_PI
#else
#define  PI	3.14159265358979323846
#endif
#undef DEG
#define	DEG	(PI/180.)


double  s_latitude = 0.66;	/* site latitude (radians north of equator) */
double  s_longitude = 2.13;	/* site longitude (radians west of Greenwich) */
double  s_meridian = 120.*DEG;	/* standard meridian (radians west) */


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
	return( 0.170 * sin( (4.*PI/373.) * (jd - 80) ) -
		0.129 * sin( (2.*PI/355.) * (jd - 8) ) +
		(12./PI) * (s_meridian - s_longitude) );
}


double
sdec(		/* solar declination angle from Julian date */
	int  jd
)
{
	return( 0.4093 * sin( (2.*PI/365.) * (jd - 81) ) );
}


double
salt(	/* solar altitude from solar declination and solar time */
	double sd,
	double st
)
{
	return( asin( sin(s_latitude) * sin(sd) -
			cos(s_latitude) * cos(sd) * cos(st*(PI/12.)) ) );
}


double
sazi(	/* solar azimuth from solar declination and solar time */
	double sd,
	double st
)
{
	return( -atan2( cos(sd)*sin(st*(PI/12.)),
 			-cos(s_latitude)*sin(sd) -
 			sin(s_latitude)*cos(sd)*cos(st*(PI/12.)) ) );
}


/****************** More accurate Michalsky algorithm ****************/


/* circle normalization */
static double
norm_cir(double r, const double p)
{
	while (r < 0) r += p;
	while (r >= p) r -= p;
	return(r);
}


/* Almanac Julian date relative to noon UT on Jan 1, 2000 (fractional days) */
double
mjdate(int year, int month, int day, double hour)
{
	int	jd = jdate(month, day);
	jd += (month > 2) & !(year&3);
	jd += (year - 1949)*365 + (year - 1949)/4;
	hour += s_meridian*(12./PI);
	return(jd + hour*(1./24.) + (2432916.5-2451545.));
}


/* Solar declination (and solar time) from Almanac Julian date (fractional) */
double
msdec(double mjd, double *stp)
{					/* Ecliptic coordinates (radians) */
	double	L = norm_cir(280.460*DEG + 0.9856474*DEG*mjd, 2.*PI);
	double	g = norm_cir(357.528*DEG + 0.9856003*DEG*mjd, 2.*PI);
	double	l = L + 1.915*DEG*sin(g) + 0.020*DEG*sin(2.*g);
	double	ep = 23.439*DEG - 4e-7*DEG*mjd;
	double	sin_l = sin(l);

	if (stp) {			/* solar time requested, also? */
		double	ra = atan2(sin_l*cos(ep), cos(l));
		double	utime = 24.*(mjd - floor(mjd)) + 12.;
		double	gmst = 6.697375 + 0.0657098242*mjd + utime;
		double	lmst = gmst - s_longitude*(12./PI);

		*stp = norm_cir(lmst - ra*(12./PI) + 12., 24.);
	}
	return(asin(sin(ep)*sin_l));	/* return solar declination angle */
}
