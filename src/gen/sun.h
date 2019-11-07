/* RCSid $Id: sun.h,v 2.3 2019/11/07 23:15:07 greg Exp $ */
/*
 * Header file for solar position calculations
 */

#ifndef _RAD_SUN_H_
#define _RAD_SUN_H_

#ifdef __cplusplus
extern "C" {
#endif
					/* sun calculation constants */
extern double  s_latitude;
extern double  s_longitude;
extern double  s_meridian;

extern int jdate(int month, int day);
extern double stadj(int  jd);
extern double sdec(int  jd);
extern double salt(double sd, double st);
extern double sazi(double sd, double st);

extern double mjdate(int year, int month, int day, double hour);
extern double msdec(double mjd, double *stp);

#ifdef __cplusplus
}
#endif

#endif /* _RAD_SUN_H_ */
