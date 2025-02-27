#ifndef lint
static const char RCSid[] = "$Id: loadEPW.c,v 2.2 2025/02/27 19:00:00 greg Exp $";
#endif
/*
 * Load an EPW (or WEA) file, one data point at a time
 *
 *	G. Ward, Feb 2025
 */

#include "rtio.h"
#include <stdlib.h>
#include <ctype.h>
#include "loadEPW.h"

const char	WDname[7][10] = {
	"Sunday", "Monday", "Tuesday", "Wednesday",
	"Thursday", "Friday", "Saturday"
};

const char	MOname[12][10] = {
	"January", "February", "March", "April", "June", "July",
	"August", "September", "October", "November", "December"
};

#define EPWDATE_INIT	{.0f,0,0,0}

static const struct DataPeriod	EPWperiodInit = {
	"DEFAULT", 1, Sunday, {.0f, 1, 0, 0}, {24.f, 31, 11, 0}
};

const EPWrecord	EPWrecInit = {
	EPWDATE_INIT, "", 99.9f, 99.9f, 999.f, 999999.f, 9999.f, 9999.f, 9999.f,
	9999.f, 9999.f, 9999.f, 999999.f, 999999.f, 999999.f, 9999.f, 999.f, 999.f,
	99.f, 99.f, 9999e3f, 99999.f, 999e-3f, 999e-3f, 999e-2f, 99,
	.0f, 1.5e-3f, .0f
};

#define EPWdateInit	(EPWrecInit.date)

static const short  mo_da[13] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

/* Get day-of-week from name (may be abbreviated) */
static int
get_day(const char dnm[])
{
	int	n = strlen(dnm);
	int	i;

	if (n < 3)
		return(-1);

	for (i = 0; i < 7; i++)
		if (!strncasecmp(dnm, WDname[i], n))
			return(i);
	return(-1);
}

/* Hour difference calculator dt1 - dt2 (crude) */
static double
hour_diff(const EPWdate *dt1, const EPWdate *dt2)
{
	double	hrdiff = dt1->hour - dt2->hour;

	hrdiff += 24.*(mo_da[dt1->month] + dt1->day - mo_da[dt2->month] - dt2->day);

	if ((dt1->year > 0) & (dt2->year > 0))
		hrdiff += 365.*24.*(dt1->year - dt2->year);

	return(hrdiff);
}

/* Check that a date/time is legal */
static int
date_OK(const EPWdate *dt)
{
	if ((0 > dt->month) | (dt->month > 11) |
			(0 > dt->hour) | (dt->hour > 24) | (dt->day < 1))
		return(0);

	return(dt->day <= (dt->month==1 ? 29 :
				mo_da[dt->month+1] - mo_da[dt->month]));
}

/* Parse date (less hour) */
static int
get_date(EPWdate *dt, const char date[])
{
	int	mo, da, yr;
	int	nf;

	*dt = EPWdateInit;
	switch (sscanf(date, "%d/%d/%d", &mo, &da, &yr)) {
	case 1:
		da = mo;	/* jdate */
		if (da <= 0)
			return(0);
		for (mo = 0; mo_da[mo+1] < da; mo++)
			if (mo >= 11)
				return(0);
		dt->day = da - mo_da[dt->month = mo];
		return(1);
	case 3:			/* month/day/year */
		dt->year = yr;
		/* fall through */
	case 2:			/* month/day */
		dt->month = mo-1;
		dt->day = da;
		return(date_OK(dt));
	}
	return(0);
}

/* Check that site data is sensible */
static int
bad_location(const struct Location *lp)
{
	if ((-90. > lp->latitude) | (lp->latitude > 90.)) {
		fputs("Bad latitude in ", stderr);
		return(1);
	}
	if ((-180. > lp->longitude) | (lp->longitude > 180.)) {
		fputs("Bad longitude in ", stderr);
		return(2);
	}
	if ((-12. > lp->timezone) | (lp->timezone > 12.)) {
		fputs("Bad time zone in ", stderr);
		return(3);
	}
	if ((-1000. > lp->elevation) | (lp->elevation > 9999.9)) {
		fputs("Bad elevation in ", stderr);
		return(4);
	}
	return(0);
}

/* Copy word to next comma in string */
static char *
tocomma(char *dst, char *src)
{
	while (*src && *src != ',')
		*dst++ = *src++;
	src += (*src == ',');
	*dst = '\0';
	return(src);
}

/* Use fscanf() to load next date and time from input */
static int
scan_date(EPWheader *epw)
{
	int	hour;
	int	minute;

	++epw->lino;
	if (epw->isWEA) {		/* simpler for WEA input */
		epw->dtpos = EPWdateInit;
		if (fscanf(epw->fp, "%hd %hd %f", &epw->dtpos.month, &epw->dtpos.day,
				&epw->dtpos.hour) != 3)
			goto scanerr;
		epw->dtpos.month--;
	} else {
					/* EPW input line */
		if (fscanf(epw->fp, "%hd,%hd,%hd,%d,%d,", &epw->dtpos.year,
				&epw->dtpos.month, &epw->dtpos.day, &hour, &minute) != 5)
			goto scanerr;
		epw->dtpos.hour = hour-1;
		if (epw->period[0].nperhour == 1)
			epw->dtpos.hour += .5;
		else
			epw->dtpos.hour += minute*(1./60.);
		epw->dtpos.month--;
	}
					/* check date/time */
	if (!date_OK(&epw->dtpos)) {
		fputs("Illegal date/time in ", stderr);
		return(0);
	}
	return(1);
scanerr:
	if (!feof(epw->fp))
		fputs("Unable to scan date in ", stderr);
	return(0);
}

/* Open input file and load header info. */
EPWheader *
EPWopen(const char *fname)
{
	char		linbuf[1024], wbuf[64];
	char		*cp;
	int		n;
	EPWheader	*hdr = (EPWheader *)calloc(1, sizeof(EPWheader));

	if (hdr == NULL)
		goto memerr;
	if (fname == NULL || !*fname) {
		fname = "<stdin>";
		hdr->fp = stdin;
	} else if ((hdr->fp = fopen(fname, "r")) == NULL) {
		perror(fname);
		free(hdr);
		return(NULL);
	}
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (!strncasecmp(linbuf, "place ", 6))
		hdr->isWEA = 1;
	else if (strncasecmp(linbuf, "LOCATION,", 9))
		goto badformat;
	if (hdr->isWEA) {		/* getting WEA header */
		cp = linbuf+6;
		if (sscanf(cp, "%[^_]_%[^\r\n]q",
				hdr->loc.city, hdr->loc.country) != 2)
			goto badformat;
		if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
			goto readerr;
		++hdr->lino;
		if (sscanf(linbuf, "latitude %lf", &hdr->loc.latitude) != 1)
			goto badformat;
		if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
			goto readerr;
		++hdr->lino;
		if (sscanf(linbuf, "longitude %lf", &hdr->loc.longitude) != 1)
			goto badformat;
		hdr->loc.longitude *= -1.;
		if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
			goto readerr;
		++hdr->lino;
		if (sscanf(linbuf, "time_zone %f", &hdr->loc.timezone) != 1)
			goto badformat;
		hdr->loc.timezone *= -1./15.;
		if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
			goto readerr;
		++hdr->lino;
		if (sscanf(linbuf, "site_elevation %f", &hdr->loc.elevation) != 1)
			goto badformat;
		if (bad_location(&hdr->loc))
			goto badformat;
		if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
			goto readerr;
		++hdr->lino;
		if (sscanf(linbuf, "weather_data_file_units %hd", &hdr->isWEA) != 1)
			goto badformat;
		hdr->nperiods = 1;
		hdr->period[0] = EPWperiodInit;
		hdr->dstart = (int)ftell(hdr->fp);
		hdr->lin0 = hdr->lino;
		if (!scan_date(hdr))	/* get date */
			goto badformat;
		return(hdr);
	}				/* else EPW header */
	cp = linbuf+9;
	cp = tocomma(hdr->loc.city, cp);
	if (!*cp) goto badformat;
	cp = tocomma(hdr->loc.state, cp);
	if (!*cp) goto badformat;
	cp = tocomma(hdr->loc.country, cp);
	if (!*cp) goto badformat;
	cp = tocomma(hdr->loc.source, cp);
	if (!*cp) goto badformat;
	cp = tocomma(hdr->loc.wmo, cp);
	if (sscanf(cp, "%lf,%lf,%f,%f", &hdr->loc.latitude, &hdr->loc.longitude,
			&hdr->loc.timezone, &hdr->loc.elevation) != 4)
		goto badformat;
	if (bad_location(&hdr->loc))
		goto badformat;
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "DESIGN CONDITIONS,", 18))
		goto badformat;
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "TYPICAL/EXTREME PERIODS,", 24))
		goto badformat;
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "GROUND TEMPERATURES,", 20))
		goto badformat;
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "HOLIDAYS/DAYLIGHT SAVINGS,", 26))
		goto badformat;
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "COMMENTS 1,", 11))
		goto badformat;
	n = strlen(linbuf+11);
	if (n > 1) {
		hdr->comments1 = (char *)malloc(n);
		if (hdr->comments1 == NULL)
			goto memerr;
		memcpy(hdr->comments1, linbuf+11, n);
		hdr->comments1[n] = '\0';
	}
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "COMMENTS 2,", 11))
		goto badformat;
	n = strlen(linbuf+11);
	if (n > 1) {
		hdr->comments2 = (char *)malloc(n);
		if (hdr->comments2 == NULL)
			goto memerr;
		memcpy(hdr->comments2, linbuf+11, n);
		hdr->comments2[n] = '\0';
	}
	if (!fgets(linbuf, sizeof(linbuf), hdr->fp))
		goto readerr;
	++hdr->lino;
	if (strncasecmp(linbuf, "DATA PERIODS,", 13))
		goto badformat;
	hdr->nperiods = 1;		/* XXX one for now */
	if (atoi(linbuf+13) != 1)
		fprintf(stderr, "%s: warning - ignoring data periods after first\n",
				fname);
	cp = strchr(linbuf+13, ',');
	if (!cp++) goto badformat;
	if ((hdr->period[0].nperhour = atoi(cp)) <= 0)
		goto badformat;
	cp = strchr(cp, ',');
	if (!cp++) goto badformat;
	cp = tocomma(hdr->period[0].name, cp);
	if (!*cp) goto badformat;
	cp = tocomma(wbuf, cp);
	if ((hdr->period[0].startday = get_day(wbuf)) < 0)
		goto badformat;
	if (!*cp) goto badformat;
	cp = tocomma(wbuf, cp);
	if (!get_date(&hdr->period[0].firstdate, wbuf))
		goto badformat;
	if (!*cp) goto badformat;
	cp = tocomma(wbuf, cp);
	if (!get_date(&hdr->period[0].lastdate, wbuf))
		goto badformat;
	hdr->period[0].lastdate.hour = 24.f;
	hdr->dstart = (int)ftell(hdr->fp);
	hdr->lin0 = hdr->lino;
	if (!scan_date(hdr))		/* and get/check date */
		goto badformat;
	if ((hdr->dtpos.month != hdr->period[0].firstdate.month) |
			(hdr->dtpos.day != hdr->period[0].firstdate.day))
		fprintf(stderr, "%s: warning - starting date does not match first record\n",
				fname);
	return(hdr);			/* seems we're good! */
memerr:
	perror("malloc");
	if (hdr) {
		if (hdr->fp != stdin) fclose(hdr->fp);
		free(hdr);
	}
	return(NULL);
readerr:
	fprintf(stderr, "%s: header read error at line %d\n", fname, hdr->lino);
	EPWclose(hdr);
	return(NULL);
badformat:
	fprintf(stderr, "%s: header format error in %s file at line %d\n",
			fname, hdr->isWEA ? "WEA" : "EPW", hdr->lino);
	EPWclose(hdr);
	return(NULL);
}

/* Seek to a particular data record if we can */
int
EPWseek(EPWheader *epw, const EPWdate *dt)
{
	if (!epw | !dt || !epw->fp || !date_OK(dt))
		return(EOF);
					/* need to back up? */
	if (feof(epw->fp) || hour_diff(&epw->dtpos, dt) > 0.99) {
		if (epw->dstart < 0 || fseek(epw->fp, epw->dstart, SEEK_SET) < 0) {
			if (epw->fp == stdin)
				fputs("Cannot seek on standard input\n", stderr);
			else
				fprintf(stderr, "Seek error on %s input file\n",
						epw->isWEA ? "WEA" : "EPW");
			return(EOF);
		}
		epw->lino = epw->lin0;
		scan_date(epw);		/* assume this much works! */
	}
					/* step through earlier records */
	while (hour_diff(&epw->dtpos, dt) < -.99) {
		int	c;
		while ((c = getc(epw->fp)) != EOF)
			if (c == '\n')	/* skip record data */
				break;
		if (scan_date(epw))	/* load next date/time */
			continue;
		if (feof(epw->fp))
			fputs("Requested date/time not found in ", stderr);
		fputs(epw->isWEA ? "WEA" : "EPW", stderr);
		fputs(" input\n", stderr);
		return(0);
	}
	return(1);
}

/* Read the next data record from input */
int
EPWread(EPWheader *epw, EPWrecord *rp)
{
	char		linbuf[1024];
	char		*cp;

	if (!epw | !rp || feof(epw->fp))
		return(EOF);

	*rp = EPWrecInit;		/* initialize to defaults */
	rp->date = epw->dtpos;		/* date already read in */
	if (!fgets(linbuf, sizeof(linbuf), epw->fp)) {
		fputs("Unexpected EOF on input\n", stderr);
		return(0);
	}
	switch (epw->isWEA) {	/* check for WEA input */
	case WEAnot:		/* nope */
		break;
	case WEAradnorm:	/* radiometric w/ direct normal */
		if (sscanf(linbuf, "%f %f", &rp->dirirrad,
				&rp->horizdiffirrad) != 2)
			goto badformat;
		break;
	case WEAradhoriz:	/* radiometric w/ direct horizontal */
		if (sscanf(linbuf, "%f %f", &rp->globhorizirrad,
				&rp->horizdiffirrad) != 2)
			goto badformat;
		rp->globhorizirrad += rp->horizdiffirrad;
		break;
	case WEAphotnorm:	/* photometric w/ direct normal */
		if (sscanf(linbuf, "%f %f", &rp->dirillum,
				&rp->diffillum) != 2)
			goto badformat;
		break;
	default:
		fputs("Illegal WEA data type\n", stderr);
		return(0);
	}
	if (epw->isWEA) {
		if (sscanf(linbuf, "%*f %*f %f %f",
				&rp->optdepth, &rp->skycover) == 2) {
			rp->optdepth *= 1e-3;
		}
		if (!scan_date(epw) && !feof(epw->fp))
			goto badformat;
		return(1);
	}
				/* continue here if EPW file */
	cp = tocomma(rp->uncert, linbuf);
	if (!*cp) goto badformat;
	if (*cp != ',') {
		rp->dbtemp = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->dptemp = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->humidity = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->atmospressure = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->exthorizirrad = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->extdirirrad = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->horizirrad_ir = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->globhorizirrad = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->dirirrad = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->horizdiffirrad = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->globhorizillum = atof(cp);
		if (rp->globhorizillum >= 999900.f)
			rp->globhorizillum = EPWrecInit.globhorizillum;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->diffillum = atof(cp);
		if (rp->diffillum >= 999900.f)
			rp->diffillum = EPWrecInit.diffillum;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->dirillum = atof(cp);
		if (rp->dirillum >= 999900.f)
			rp->dirillum = EPWrecInit.dirillum;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->zenlum = atof(cp);
		if (rp->zenlum >= 9999.f)
			rp->zenlum = EPWrecInit.zenlum;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->windirection = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->windspeed = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->skycover = atof(cp) * 0.1;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->opskycover = atof(cp) * 0.1;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->visibility = atof(cp) * 1000.;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->ceilheight = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	cp = strchr(cp, ',');
	if (!cp++) goto badformat;
	cp = strchr(cp, ',');
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->precip = atof(cp) * 1e-3;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->optdepth = atof(cp) * 1e-3;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->snowdepth = atof(cp) * 1e-2;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->nosnowdays = atoi(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->albedo = atof(cp);
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if (*cp != ',') {
		rp->liqpdepth = atof(cp) * 1e-3;
		cp = strchr(cp, ',');
	}
	if (!cp++) goto badformat;
	if ((*cp != ',') & (*cp != '\n'))
		rp->liqhours = atof(cp);
	if (scan_date(epw) || feof(epw->fp))
		return(1);	/* normal return (even if next is EOF) */
badformat:
	fprintf(stderr, "%s file, format error at line %d\n",
			epw->isWEA ? "WEA" : "EPW", epw->lino);
	return(0);
}

/* Close input and free header data */
void
EPWclose(EPWheader *epw)
{
	if (!epw) return;
	if (epw->fp != stdin) fclose(epw->fp);
	if (epw->comments1) free(epw->comments1);
	if (epw->comments2) free(epw->comments2);
	free(epw);
	return;
}

#ifdef TEST_MAIN
int
main(int argc, char *argv[])
{
	EPWheader	*epw = EPWopen(argv[1]);
	EPWrecord	erec;
	int		rval;

	if (!epw) return(1);

	printf("Weather input file type is %d\n", epw->isWEA);
	printf("Location is %s, %s\n", epw->loc.city, epw->loc.country);
	printf("Time zone is %.2f hours from GMT\n", epw->loc.timezone);
	printf("Latitude, longitude is (%.8f,%.8f) degrees (N,E)\n",
			epw->loc.latitude, epw->loc.longitude);
	printf("Elevation is %.0f meters\n", epw->loc.elevation);
	printf("'%s' starts %d/%d on a %s and ends %d/%d\n", epw->period[0].name,
			epw->period[0].firstdate.month+1, epw->period[0].firstdate.day,
			WDname[epw->period[0].startday],
			epw->period[0].lastdate.month+1, epw->period[0].lastdate.day);
	while ((rval = EPWread(epw, &erec)) > 0) {
		printf("Weather record for %d/%d/%d @ %.2f hours:\n", erec.date.month+1,
				erec.date.day, erec.date.year, erec.date.hour);
		if (EPWisset(&erec,dbtemp))
			printf("\tDry bulb temp: %.1f °C\n", erec.dbtemp);
		if (EPWisset(&erec,dptemp))
			printf("\tDew point temp: %.1f °C\n", erec.dptemp);
		if (EPWisset(&erec,humidity))
			printf("\tHumidity: %.1f%%\n", erec.humidity);
		if (EPWisset(&erec,atmospressure))
			printf("\tAtmospheric pressure: %.0f Pa\n", erec.atmospressure);
		if (EPWisset(&erec,exthorizirrad))
			printf("\tExtraterrestrial horiz irrad: %.1f w/m2\n", erec.exthorizirrad);
		if (EPWisset(&erec,horizirrad_ir))
			printf("\tInfrared direct normal irrad: %.1f w/m2\n", erec.horizirrad_ir);
		if (EPWisset(&erec,globhorizirrad))
			printf("\tGlobal horiz irrad: %.1f w/m2\n", erec.globhorizirrad);
		if (EPWisset(&erec,dirirrad))
			printf("\tDirect normal irrad: %.1f w/m2\n", erec.dirirrad);
		if (EPWisset(&erec,horizdiffirrad))
			printf("\tHorizontal diffuse irrad: %.1f w/m2\n", erec.horizdiffirrad);
		if (EPWisset(&erec,globhorizillum))
			printf("\tGlobal horizontal illuminance: %.0f lux\n", erec.globhorizillum);
		if (EPWisset(&erec,diffillum))
			printf("\tDiffuse horizontal illuminance: %.0f lux\n", erec.diffillum);
		if (EPWisset(&erec,dirillum))
			printf("\tDirect normal illuminance: %.0f lux\n", erec.dirillum);
		if (EPWisset(&erec,zenlum))
			printf("\tZenith luminance: %.0f nits\n", erec.zenlum);
		if (EPWisset(&erec,windirection))
			printf("\tWind direction: %.1f °E of S\n", erec.windirection);
		if (EPWisset(&erec,windspeed))
			printf("\tWind speed: %.2f m/s\n", erec.windspeed);
		if (EPWisset(&erec,skycover))
			printf("\tSky cover fraction: %.3f\n", erec.skycover);
		if (EPWisset(&erec,opskycover))
			printf("\tOpaque sky cover fraction: %.3f\n", erec.opskycover);
		if (EPWisset(&erec,visibility))
			printf("\tVisibility: %.0f meters\n", erec.visibility);
		if (EPWisset(&erec,ceilheight))
			printf("\tCloud ceiling height: %.0f meters\n", erec.ceilheight);
		if (EPWisset(&erec,precip))
			printf("\tPrecipitation: %f meters\n", erec.precip);
		if (EPWisset(&erec,optdepth))
			printf("\tAerosol optical depth: %.4f\n", erec.optdepth);
		if (EPWisset(&erec,snowdepth))
			printf("\tSnow depth: %.2f meters\n", erec.snowdepth);
		if (EPWisset(&erec,nosnowdays))
			printf("\tDays since last snow: %d\n", erec.nosnowdays);
		if (EPWisset(&erec,albedo))
			printf("\tAlbedo: %.4f\n", erec.albedo);
		if (EPWisset(&erec,liqpdepth))
			printf("\tLiquid precipitation depth: %.6f meters\n", erec.liqpdepth);
		if (EPWisset(&erec,liqhours))
			printf("\tLiquid precipitation time: %.1f hours\n", erec.liqhours);
	}
	EPWclose(epw);
	return(!rval);
}
#endif	/* TEST_MAIN */
