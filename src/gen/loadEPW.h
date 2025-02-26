/* RCSid $Id: loadEPW.h,v 2.1 2025/02/26 20:39:28 greg Exp $ */
/*
 * Header for EPW (or WEA) file loader
 *
 *	G. Ward, Feb 2025
 */

#ifndef _RAD_LOADEPW_H_
#define _RAD_LOADEPW_H_

#ifdef __cplusplus
extern "C" {
#endif

enum {Sunday, Monday, Tuesday, Wednesday,
	Thursday, Friday, Saturday};

extern const char	WDname[7][10];

enum {January, February, March, April, June, July,
	August, September, October, November, December};

extern const char	MOname[12][10];

enum {WEAnot=0, WEAradnorm=1, WEAradhoriz=2, WEAphotnorm=3};

/* EPW date and time */
typedef struct {
	float	hour;		/* decimal hour (0-24) */
	short	day;		/* day of the month (1-31) */
	short	month;		/* month of the year (0-11) */
	short	year;		/* year (or 0) */
} EPWdate;

/* EPW file header info. */
typedef struct {
	FILE	*fp;			/* open file pointer */
	int	lino;			/* current line in file */
	int	dstart;			/* starting data position in file */
	short	lin0;			/* # lines in header */
	short	isWEA;			/* input is WEA data type? */
	EPWdate	dtpos;			/* date and time of current record */
	struct Location {
		char	city[64];		/* location city */
		char	state[64];		/* location state/province */
		char	country[64];		/* country */
		char	source[64];		/* data source */
		char	wmo[16];		/* WMO (EnergyPlus alpha?) */
		double	latitude;		/* site latitude (deg north) */
		double	longitude;		/* site longitude (deg east) */
		float	timezone;		/* time zone (hours from GMT) */
		float	elevation;		/* site altitude in meters */
	} loc;				/* location data */
	/* struct DesConditions		XXX unsupported */
	/* struct VarPeriods		XXX unsupported */
	/* struct GroundTemps		XXX unsupported */
	/* struct SpecialDays		XXX unsupported */
	char	*comments1;		/* first set of additional comments */
	char	*comments2;		/* second set of comments */
	int	nperiods;		/* # of data periods */
	struct DataPeriod {
		char	name[64];	/* data period name */
		short	nperhour;	/* # records/hour */
		short	startday;	/* starting day of the week */
		EPWdate	firstdate;	/* starting date */
		EPWdate	lastdate;	/* ending date */
	} period[1];		/* periods covered in EPW file */
} EPWheader;

/* EPW record holder (NOTE: all lengths in meters, contrary to spec!) */
typedef struct {
	EPWdate		date;		/* measurement date and time */
	char		uncert[64];	/* uncertainty flags */
	float		dbtemp;		/* dry bulb temperature (deg-C) */
	float		dptemp;		/* dew point temperature (deg-C) */
	float		humidity;	/* relative humidity (0-110 %) */
	float		atmospressure;	/* atmospheric pressure (31000 - 120000) */
	float		exthorizirrad;	/* extraterrestrial horizontal irradiance */
	float		extdirirrad;	/* extraterrestrial direct normal irradiance */
	float		horizirrad_ir;	/* infrared horizontal normal irradiance */
	float		globhorizirrad;	/* global horizontal irradiance */
	float		dirirrad;	/* direct normal irradiance */
	float		horizdiffirrad;	/* horizontal diffuse irradiance */
	float		globhorizillum;	/* global horizontal illuminance (lux) */
	float		diffillum;	/* diffuse horizontal illuminance (lux) */
	float		dirillum;	/* direct normal illuminance (lux) */
	float		zenlum;		/* zenith luminance (cd/m2) */
	float		windirection;	/* wind direction (deg E of N) */
	float		windspeed;	/* wind speed in m/s */
	float		skycover;	/* total sky coverage (fraction) */
	float		opskycover;	/* opaque sky coverage (fraction) */
	float		visibility;	/* visibility (meters!) */
	float		ceilheight;	/* cloud ceiling height (meters) */
	/* char		weather[32];	XXX unsupported */
	float		precip;		/* precipital water in meters(!) */
	float		optdepth;	/* aerosol optical depth (fraction!) */
	float		snowdepth;	/* snow depth (meters!) */
	int		nosnowdays;	/* # days since last snow */
	float		albedo;		/* ratio of reflected to global horiz. irrad */
	float		liqpdepth;	/* liquid precipitation depth (meters!) */
	float		liqhours;	/* liquid precipitation quantity (hours) */
} EPWrecord;

extern const EPWrecord	EPWrecInit;

#define EPWisset(rp,field)	((rp)->field != EPWrecInit.field)

				/* specify NULL for stdin */
extern EPWheader *	EPWopen(const char *fname);

extern int		EPWseek(EPWheader *epw, const EPWdate *dt);

extern int		EPWread(EPWheader *epw, EPWrecord *rp);

extern void		EPWclose(EPWheader *epw);

#ifdef __cplusplus
}
#endif
#endif		/* ! _RAD_LOADEPW_H_ */
