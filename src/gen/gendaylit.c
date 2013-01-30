#ifndef lint
static const char RCSid[] = "$Id: gendaylit.c,v 2.9 2013/01/30 01:02:42 greg Exp $";
#endif
/*      Copyright (c) 1994,2006 *Fraunhofer Institut for Solar Energy Systems
 *        			Heidenhofstr. 2, D-79110 Freiburg, Germany
 *      			*Agence de l'Environnement et de la Maitrise de l'Energie
 *      			Centre de Valbonne, 500 route des Lucioles, 06565 Sophia Antipolis Cedex, France
 * 				*BOUYGUES 
 *				1 Avenue Eugene Freyssinet, Saint-Quentin-Yvelines, France
*
*      24.1.2006                some adjustments for cygwin compilation, inclusion of RADIANCE3.7 libraries, by J. Wienold
*		2011/10/08	chr@riap.de:
*						- integrated coeff_perez.dat and defangles.dat
*						- avoid some segfaults caused by out of range parameters and
*						- numerically dangerous range checks 
*/

/*
 *  gendaylit.c 	program to generate the angular distribution of the daylight.
 *			Our zenith is along the Z-axis, the X-axis
 *			points east, and the Y-axis points north.
*/

#include  <stdio.h>
#include  <string.h>
#include  <math.h>
#include  <stdlib.h>

#include  "color.h"
#include  "paths.h"

#define  DOT(v1,v2)	(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

double  normsc();

/*static	char *rcsid="$Header: /home/cvsd/radiance/ray/src/gen/gendaylit.c,v 2.9 2013/01/30 01:02:42 greg Exp $";*/

float coeff_perez[] = {
	1.3525,-0.2576,-0.2690,-1.4366,-0.7670,0.0007,1.2734,-0.1233,2.8000,0.6004,1.2375,1.000,1.8734,0.6297,0.9738,0.2809,0.0356,-0.1246,-0.5718,0.9938,
	-1.2219,-0.7730,1.4148,1.1016,-0.2054,0.0367,-3.9128,0.9156,6.9750,0.1774,6.4477,-0.1239,-1.5798,-0.5081,-1.7812,0.1080,0.2624,0.0672,-0.2190,-0.4285,
	-1.1000,-0.2515,0.8952,0.0156,0.2782,-0.1812,-4.5000,1.1766,24.7219,-13.0812,-37.7000,34.8438,-5.0000,1.5218,3.9229,-2.6204,-0.0156,0.1597,0.4199,-0.5562,
	-0.5484,-0.6654,-0.2672,0.7117,0.7234,-0.6219,-5.6812,2.6297,33.3389,-18.3000,-62.2500,52.0781,-3.5000,0.0016,1.1477,0.1062,0.4659,-0.3296,-0.0876,-0.0329,
	-0.6000,-0.3566,-2.5000,2.3250,0.2937,0.0496,-5.6812,1.8415,21.0000,-4.7656,-21.5906,7.2492,-3.5000,-0.1554,1.4062,0.3988,0.0032,0.0766,-0.0656,-0.1294,
	-1.0156,-0.3670,1.0078,1.4051,0.2875,-0.5328,-3.8500,3.3750,14.0000,-0.9999,-7.1406,7.5469,-3.4000,-0.1078,-1.0750,1.5702,-0.0672,0.4016,0.3017,-0.4844,
	-1.0000,0.0211,0.5025,-0.5119,-0.3000,0.1922,0.7023,-1.6317,19.0000,-5.0000,1.2438,-1.9094,-4.0000,0.0250,0.3844,0.2656,1.0468,-0.3788,-2.4517,1.4656,
	-1.0500,0.0289,0.4260,0.3590,-0.3250,0.1156,0.7781,0.0025,31.0625,-14.5000,-46.1148,55.3750,-7.2312,0.4050,13.3500,0.6234,1.5000,-0.6426,1.8564,0.5636};


float defangle_theta[] = {84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 12, 12, 12, 12, 12, 12, 0};

float defangle_phi[] = {0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240, 252, 264, 276, 288, 300, 312, 324, 336, 348, 0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240, 252, 264, 276, 288, 300, 312, 324, 336, 348, 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 0, 60, 120, 180, 240, 300, 0};



/* Perez sky parametrization : epsilon and delta calculations from the direct and diffuse irradiances */
double sky_brightness();
double sky_clearness();

/* calculation of the direct and diffuse components from the Perez parametrization */
double	diffuse_irradiance_from_sky_brightness();
double 	direct_irradiance_from_sky_clearness();


/* Perez global horizontal, diffuse horizontal and direct normal luminous efficacy models : input w(cm)=2cm, solar zenith angle(degrees); output efficacy(lm/W) */
double 	glob_h_effi_PEREZ();
double 	glob_h_diffuse_effi_PEREZ();
double 	direct_n_effi_PEREZ();
/*likelihood check of the epsilon, delta, direct and diffuse components*/
void 	check_parametrization();
void 	check_irradiances();
void 	check_illuminances();
void 	illu_to_irra_index();
void	print_error_sky();


/* Perez sky luminance model */
double  calc_rel_lum_perez(double dzeta,double gamma,double Z,
		double epsilon,double Delta,float coeff_perez[]);
/* coefficients for the sky luminance perez model */
void 	coeff_lum_perez(double Z, double epsilon, double Delta, float coeff_perez[]);
double	radians(double degres);
double	degres(double radians);
void	theta_phi_to_dzeta_gamma(double theta,double phi,double *dzeta,double *gamma, double Z);
double 	integ_lv(float *lv,float *theta);

void printdefaults();
void computesky();
void printhead(int ac, char** av);
void userror(char* msg);
void printsky();

FILE * frlibopen(char* fname);

/* astronomy and geometry*/
double 	get_eccentricity();
double 	air_mass();

extern int jdate(int month, int day);
extern double stadj(int  jd);
extern double sdec(int  jd);
extern double salt(double sd, double st);
extern double sazi(double sd, double st);


/* sun calculation constants */
extern double  s_latitude;
extern double  s_longitude;
extern double  s_meridian;

const double	AU = 149597890E3;
const double 	solar_constant_e = 1367;    /* solar constant W/m^2 */
const double  	solar_constant_l = 127.5;   /* solar constant klux */

const double	half_sun_angle = 0.2665;
const double	half_direct_angle = 2.85;

const double	skyclearinf = 1.000;	/* limitations for the variation of the Perez parameters */
const double	skyclearsup = 12.1;
const double	skybriginf = 0.01;
const double	skybrigsup = 0.6;



/* required values */
int  	month, day;				/* date */
double  hour;					/* time */
int  	tsolar;					/* 0=standard, 1=solar */
double  altitude, azimuth;			/* or solar angles */



/* definition of the sky conditions through the Perez parametrization */
double	skyclearness = 0; 
double  skybrightness = 0;
double	solarradiance;	/*radiance of the sun disk and of the circumsolar area*/
double	diffuseilluminance, directilluminance, diffuseirradiance, directirradiance;
double	sunzenith, daynumber=150, atm_preci_water=2;

double sunaltitude_border = 0;
double 	diffnormalization = 0;
double dirnormalization = 0;
double 	*c_perez;

int	output=0;	/*define the unit of the output (sky luminance or radiance): visible watt=0, solar watt=1, lumen=2*/
int	input=0;	/*define the input for the calulation*/

int	suppress_warnings=0;

	/* default values */
int  cloudy = 0;				/* 1=standard, 2=uniform */
int  dosun = 1;
double  zenithbr = -1.0;
double  betaturbidity = 0.1;
double  gprefl = 0.2;
int	S_INTER=0;

	/* computed values */
double  sundir[3];
double  groundbr = 0;
double  F2;
double  solarbr = 0.0;
int	u_solar = 0;				/* -1=irradiance, 1=radiance */

char  *progname;
char  errmsg[128];


int main(int argc, char** argv)
{
	int  i;

	progname = argv[0];
	if (argc == 2 && !strcmp(argv[1], "-defaults")) {
		printdefaults();
		return 0;
	}
	if (argc < 4)
		userror("arg count");
	if (!strcmp(argv[1], "-ang")) {
		altitude = atof(argv[2]) * (M_PI/180);
		azimuth = atof(argv[3]) * (M_PI/180);
		month = 0;
	} else {
		month = atoi(argv[1]);
		if (month < 1 || month > 12)
			userror("bad month");
		day = atoi(argv[2]);
		if (day < 1 || day > 31)
			userror("bad day");
		hour = atof(argv[3]);
		if (hour < 0 || hour >= 24)
			userror("bad hour");
		tsolar = argv[3][0] == '+';
	}
	for (i = 4; i < argc; i++)
		if (argv[i][0] == '-' || argv[i][0] == '+')
			switch (argv[i][1]) {
			case 's':
				cloudy = 0;
				dosun = argv[i][0] == '+';
				break;
			case 'r':
			case 'R':
				u_solar = argv[i][1] == 'R' ? -1 : 1;
				solarbr = atof(argv[++i]);
				break;
			case 'c':
				cloudy = argv[i][0] == '+' ? 2 : 1;
				dosun = 0;
				break;
			case 't':
				betaturbidity = atof(argv[++i]);
				break;
			case 'w':
				suppress_warnings = 1;
				break;			
			case 'b':
				zenithbr = atof(argv[++i]);
				break;
			case 'g':
				gprefl = atof(argv[++i]);
				break;
			case 'a':
				s_latitude = atof(argv[++i]) * (M_PI/180);
				break;
			case 'o':
				s_longitude = atof(argv[++i]) * (M_PI/180);
				break;
			case 'm':
				s_meridian = atof(argv[++i]) * (M_PI/180);
				break;
			
			case 'O':
				output = atof(argv[++i]);	/*define the unit of the output of the program : 
								sky and sun luminance/radiance (0==W visible, 1==W solar radiation, 2==lm) */
				break;
				
			case 'P':
				input = 0;				/* Perez parameters: epsilon, delta */
				skyclearness = atof(argv[++i]);
				skybrightness = atof(argv[++i]);
				break;

			case 'W':					/* direct normal Irradiance [W/m^2] */
				input = 1;				/* diffuse horizontal Irrad. [W/m^2] */
				directirradiance = atof(argv[++i]);
				diffuseirradiance = atof(argv[++i]);
				break;
				
			case 'L':					/* direct normal Illuminance [Lux] */
				input = 2;				/* diffuse horizontal Ill. [Lux] */
				directilluminance = atof(argv[++i]);
				diffuseilluminance = atof(argv[++i]);
				break;
			
			case 'G':					/* direct horizontal Irradiance [W/m^2] */
				input = 3;				/* diffuse horizontal Irrad. [W/m^2] */
				directirradiance = atof(argv[++i]);
				diffuseirradiance = atof(argv[++i]);
				break;
			
		 	case 'l':
				sunaltitude_border = atof(argv[++i]);
				break;
			
			
			default:
				sprintf(errmsg, "unknown option: %s", argv[i]);
				userror(errmsg);
			}
		else
			userror("bad option");

	if (fabs(s_meridian-s_longitude) > 30*M_PI/180)
		fprintf(stderr,
		    "%s: warning: %.1f hours btwn. standard meridian and longitude\n",
		    progname, (s_longitude-s_meridian)*12/M_PI);


	/* allocation dynamique de memoire pour les pointeurs */
	if ( (c_perez = calloc(5, sizeof(double))) == NULL )
	{
		fprintf(stderr,"Out of memory error in function main !");
		return 1;
	}

	printhead(argc, argv);

	computesky();
	printsky();

	return 0;
}


void computesky()			/* compute sky parameters */
{

	/* new variables */
	int 	j;
	float 	*lv_mod;  /* 145 luminance values*/
	  /* 145 directions for the calculation of the normalization coefficient, coefficient Perez model */
	float 	*theta_o, *phi_o;
	double 	dzeta, gamma;
	double	normfactor;



	/* compute solar direction */

	if (month) {			/* from date and time */
		int  jd;
		double  sd, st;

		jd = jdate(month, day);		/* Julian date */
		sd = sdec(jd);			/* solar declination */
		if (tsolar)			/* solar time */
			st = hour;
		else
			st = hour + stadj(jd);
		altitude = salt(sd, st);
		azimuth = sazi(sd, st);
		
		daynumber = (double)jdate(month, day);

	}
	
	
	
	
	
	/* if loop for the -l option. 01/2013 Sprenger  */
	
	if (altitude*180/M_PI < sunaltitude_border) {
	
	if (suppress_warnings==0)
	    fprintf(stderr, "Warning: sun altitude (%.3f degrees) below the border (%.3f degrees)\n",altitude*180/M_PI,sunaltitude_border);
	print_error_sky();
	exit(0);
	}
	
	
			
	
			
	if (!cloudy && altitude > 87.*M_PI/180.) {
		
		if (suppress_warnings==0) { 
		    fprintf(stderr,
		    "%s: warning - sun too close to zenith, reducing altitude to 87 degrees\n",
		    progname);
		}
		altitude = 87.*M_PI/180.;
	}
	
	sundir[0] = -sin(azimuth)*cos(altitude);
	sundir[1] = -cos(azimuth)*cos(altitude);
	sundir[2] = sin(altitude);
	
		
	/* calculation for the new functions */
	sunzenith = 90 - altitude*180/M_PI;
	
	

        /* compute the inputs for the calculation of the light distribution over the sky*/
	if (input==0)
		{
		check_parametrization();
		diffuseirradiance = diffuse_irradiance_from_sky_brightness(); /*diffuse horizontal irradiance*/
		directirradiance = direct_irradiance_from_sky_clearness();
		check_irradiances();
		
		if (output==0 || output==2)
			{
			diffuseilluminance = diffuseirradiance*glob_h_diffuse_effi_PEREZ();/*diffuse horizontal illuminance*/
			directilluminance = directirradiance*direct_n_effi_PEREZ();
			check_illuminances();
			}
		}
	

	else if (input==1)
		{
		check_irradiances();
		skybrightness = sky_brightness();
		skyclearness =  sky_clearness();
		check_parametrization();

		if (output==0 || output==2)
			{
			diffuseilluminance = diffuseirradiance*glob_h_diffuse_effi_PEREZ();/*diffuse horizontal illuminance*/
			directilluminance = directirradiance*direct_n_effi_PEREZ();
			check_illuminances();
			}

		}
			
	
	else if (input==2)
		{		
		check_illuminances();
		illu_to_irra_index();
		check_parametrization();
		}
		

	else if (input==3)
		{
			if (altitude<=0)
			{
				if (suppress_warnings==0)
				     fprintf(stderr, "Warning: solar zenith angle larger than 90 degrees; using zero irradiance to proceed\n");
				directirradiance = 0;
				diffuseirradiance = 0;
			} else {
				directirradiance=directirradiance/sin(altitude);
			}
		check_irradiances();
		skybrightness = sky_brightness();
		skyclearness =  sky_clearness();
		check_parametrization();

		if (output==0 || output==2)
			{
			diffuseilluminance = diffuseirradiance*glob_h_diffuse_effi_PEREZ();/*diffuse horizontal illuminance*/
			directilluminance = directirradiance*direct_n_effi_PEREZ();
			check_illuminances();
			}

		}

	
	else	{fprintf(stderr,"error in giving the input arguments"); exit(1);}


	
        /* normalization factor for the relative sky luminance distribution, diffuse part*/

	if ( (lv_mod = malloc(145*sizeof(float))) == NULL)
	{
		fprintf(stderr,"Out of memory in function main");
		exit(1);
	}

	/* read the angles */
	theta_o = defangle_theta;
	phi_o = defangle_phi;

        /* parameters for the perez model */
	coeff_lum_perez(radians(sunzenith), skyclearness, skybrightness, coeff_perez);

        /*calculation of the modelled luminance */
	for (j=0;j<145;j++)
	{
		theta_phi_to_dzeta_gamma(radians(*(theta_o+j)),radians(*(phi_o+j)),&dzeta,&gamma,radians(sunzenith));
		*(lv_mod+j) = calc_rel_lum_perez(dzeta,gamma,radians(sunzenith),skyclearness,skybrightness,coeff_perez);
		// printf("theta, phi, lv_mod %f\t %f\t %f\n", *(theta_o+j),*(phi_o+j),*(lv_mod+j));
	}

	/* integration of luminance for the normalization factor, diffuse part of the sky*/
	diffnormalization = integ_lv(lv_mod, theta_o);
	/*printf("perez integration %lf\n", diffnormalization);*/
	
	


        /*normalization coefficient in lumen or in watt*/
	if (output==0)
		{
		diffnormalization = diffuseilluminance/diffnormalization/WHTEFFICACY;
		}
	else if (output==1)
		{
		diffnormalization = diffuseirradiance/diffnormalization;
		}
	else if (output==2)
		{
		diffnormalization = diffuseilluminance/diffnormalization;
		}

	else	{fprintf(stderr,"Wrong output specification.\n"); exit(1);}




        /* calculation for the solar source */	
	if (output==0) 
		solarradiance = directilluminance/(2*M_PI*(1-cos(half_sun_angle*M_PI/180)))/WHTEFFICACY;
		
	else if (output==1)
		solarradiance = directirradiance/(2*M_PI*(1-cos(half_sun_angle*M_PI/180)));
	
	else
		solarradiance = directilluminance/(2*M_PI*(1-cos(half_sun_angle*M_PI/180)));
	
			


/* Compute the ground radiance */
zenithbr=calc_rel_lum_perez(0.0,radians(sunzenith),radians(sunzenith),skyclearness,skybrightness,coeff_perez);
zenithbr*=diffnormalization;
	
if (skyclearness==1)
	normfactor = 0.777778;
		
if (skyclearness>=6)
	{		
	F2 = 0.274*(0.91 + 10.0*exp(-3.0*(M_PI/2.0-altitude)) + 0.45*sundir[2]*sundir[2]);
	normfactor = normsc()/F2/M_PI;
	}

if ( (skyclearness>1) && (skyclearness<6) )
	{
	S_INTER=1;
	F2 = (2.739 + .9891*sin(.3119+2.6*altitude)) * exp(-(M_PI/2.0-altitude)*(.4441+1.48*altitude));
	normfactor = normsc()/F2/M_PI;
	}

groundbr = zenithbr*normfactor;

if (dosun&&(skyclearness>1)) 
	groundbr += 6.8e-5/M_PI*solarradiance*sundir[2];		

groundbr *= gprefl;



return;
}




void print_error_sky()
{
	sundir[0] = -sin(azimuth)*cos(altitude);
	sundir[1] = -cos(azimuth)*cos(altitude);
	sundir[2] = sin(altitude);
	
	printf("\nvoid brightfunc skyfunc\n");
	printf("2 skybright perezlum.cal\n");
	printf("0\n");
	printf("10 0.00 0.00  0.000 0.000 0.000 0.000 0.000  %f %f %f \n", sundir[0], sundir[1], sundir[2]);
}
	


void printsky()			/* print out sky */
{
	if (dosun&&(skyclearness>1)) 
	{		
		printf("\nvoid light solar\n");
		printf("0\n0\n");
		printf("3 %.3e %.3e %.3e\n", solarradiance, solarradiance, solarradiance);
		printf("\nsolar source sun\n");
		printf("0\n0\n");
		printf("4 %f %f %f %f\n", sundir[0], sundir[1], sundir[2], 2*half_sun_angle);
	} else if (dosun) {
		printf("\nvoid light solar\n");
		printf("0\n0\n");
		printf("3 0.0 0.0 0.0\n");
		printf("\nsolar source sun\n");
		printf("0\n0\n");
		printf("4 %f %f %f %f\n", sundir[0], sundir[1], sundir[2], 2*half_sun_angle);
	}
	
	printf("\nvoid brightfunc skyfunc\n");
	printf("2 skybright perezlum.cal\n");
	printf("0\n");
	printf("10 %.3e %.3e %lf %lf %lf %lf %lf %f %f %f \n", diffnormalization, groundbr,
	    	*(c_perez+0),*(c_perez+1),*(c_perez+2),*(c_perez+3),*(c_perez+4), 
	    	sundir[0], sundir[1], sundir[2]);
}


void printdefaults()			/* print default values */
{
	printf("-g %f\t\t\t# Ground plane reflectance\n", gprefl);
	if (zenithbr > 0.0)
		printf("-b %f\t\t\t# Zenith radiance (watts/ster/m^2\n", zenithbr);
	else
		printf("-t %f\t\t\t# Atmospheric betaturbidity\n", betaturbidity);
	printf("-a %f\t\t\t# Site latitude (degrees)\n", s_latitude*(180/M_PI));
	printf("-o %f\t\t\t# Site longitude (degrees)\n", s_longitude*(180/M_PI));
	printf("-m %f\t\t\t# Standard meridian (degrees)\n", s_meridian*(180/M_PI));
}


void userror(char* msg)			/* print usage error and quit */
{
	if (msg != NULL)
		fprintf(stderr, "%s: Use error - %s\n", progname, msg);
	fprintf(stderr, "Usage: %s month day hour [-P|-W|-L|-G] direct_value diffuse_value [options]\n", progname);
	fprintf(stderr, "or:    %s -ang altitude azimuth [-P|-W|-L|-G] direct_value diffuse_value [options]\n", progname);
	fprintf(stderr, "	-P epsilon delta  (these are the Perez parameters) \n");
	fprintf(stderr, "	-W direct-normal-irradiance diffuse-horizontal-irradiance (W/m^2)\n");
	fprintf(stderr, "	-L direct-normal-illuminance diffuse-horizontal-illuminance (lux)\n");
	fprintf(stderr, "	-G direct-horizontal-irradiance diffuse-horizontal-irradiance (W/m^2)\n");
	fprintf(stderr, "	-O [0|1|2]  (0=output in W/m^2/sr visible, 1=output in W/m^2/sr solar, 2=output in candela/m^2), default is 0 \n");
	fprintf(stderr, "	gendaylit version 2.00 (2013/01/28)  \n");
	exit(1);
}



double normsc()	          /* compute normalization factor (E0*F2/L0) */
{
	static double  nfc[2][5] = {
				/* clear sky approx. */
		{2.766521, 0.547665, -0.369832, 0.009237, 0.059229},
				/* intermediate sky approx. */
		{3.5556, -2.7152, -1.3081, 1.0660, 0.60227},
	};
	register double  *nf;
	double  x, nsc;
	register int  i;
					/* polynomial approximation */
	nf = nfc[S_INTER];
	x = (altitude - M_PI/4.0)/(M_PI/4.0);
	nsc = nf[i=4];
	while (i--)
		nsc = nsc*x + nf[i];

	return(nsc);
}



void printhead(int ac, char** av)		/* print command header */
{
	putchar('#');
	while (ac--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');
}






/* Perez models */

/* Perez global horizontal luminous efficacy model */
double glob_h_effi_PEREZ()
{

	double 	value;
	double 	category_bounds[10], a[10], b[10], c[10], d[10];
	int   	category_total_number, category_number, i;


if ((skyclearness<skyclearinf || skyclearness>skyclearsup || skybrightness<skybriginf || skybrightness>skybrigsup) && suppress_warnings==0)
     fprintf(stderr, "Warning: skyclearness or skybrightness out of range in function glob_h_effi_PEREZ \n");

	/* initialize category bounds (clearness index bounds) */

	category_total_number = 8;

	category_bounds[1] = 1;
	category_bounds[2] = 1.065;
	category_bounds[3] = 1.230;
	category_bounds[4] = 1.500;
	category_bounds[5] = 1.950;
	category_bounds[6] = 2.800;
	category_bounds[7] = 4.500;
	category_bounds[8] = 6.200;
	category_bounds[9] = 12.01;


	/* initialize model coefficients */
	a[1] = 96.63;
	a[2] = 107.54;
	a[3] = 98.73;
	a[4] = 92.72;
	a[5] = 86.73;
	a[6] = 88.34;
	a[7] = 78.63;
	a[8] = 99.65;

	b[1] = -0.47;
	b[2] = 0.79;
	b[3] = 0.70;
	b[4] = 0.56;
	b[5] = 0.98;
	b[6] = 1.39;
	b[7] = 1.47;
	b[8] = 1.86;

	c[1] = 11.50;
	c[2] = 1.79;
	c[3] = 4.40;
	c[4] = 8.36;
	c[5] = 7.10;
	c[6] = 6.06;
	c[7] = 4.93;
	c[8] = -4.46;

	d[1] = -9.16;
	d[2] = -1.19;
	d[3] = -6.95;
	d[4] = -8.31;
	d[5] = -10.94;
	d[6] = -7.60;
	d[7] = -11.37;
	d[8] = -3.15;




	for (i=1; i<=category_total_number; i++)
	{
		if ( (skyclearness >= category_bounds[i]) && (skyclearness < category_bounds[i+1]) )
			category_number = i;
	}

	value = a[category_number] + b[category_number]*atm_preci_water  + 
	    c[category_number]*cos(sunzenith*M_PI/180) +  d[category_number]*log(skybrightness);

	return(value);
}


/* global horizontal diffuse efficacy model, according to PEREZ */
double glob_h_diffuse_effi_PEREZ()
{
	double 	value;
	double 	category_bounds[10], a[10], b[10], c[10], d[10];
	int   	category_total_number, category_number, i;

	
if ((skyclearness<skyclearinf || skyclearness>skyclearsup || skybrightness<skybriginf || skybrightness>skybrigsup) && suppress_warnings==0)
     fprintf(stderr, "Warning: skyclearness or skybrightness out of range in function glob_h_diffuse_effi_PEREZ \n");

/* initialize category bounds (clearness index bounds) */

	category_total_number = 8;

//XXX: category_bounds > 0.1
	category_bounds[1] = 1;
	category_bounds[2] = 1.065;
	category_bounds[3] = 1.230;
	category_bounds[4] = 1.500;
	category_bounds[5] = 1.950;
	category_bounds[6] = 2.800;
	category_bounds[7] = 4.500;
	category_bounds[8] = 6.200;
	category_bounds[9] = 12.01;


	/* initialize model coefficients */
	a[1] = 97.24;
	a[2] = 107.22;
	a[3] = 104.97;
	a[4] = 102.39;
	a[5] = 100.71;
	a[6] = 106.42;
	a[7] = 141.88;
	a[8] = 152.23;

	b[1] = -0.46;
	b[2] = 1.15;
	b[3] = 2.96;
	b[4] = 5.59;
	b[5] = 5.94;
	b[6] = 3.83;
	b[7] = 1.90;
	b[8] = 0.35;

	c[1] = 12.00;
	c[2] = 0.59;
	c[3] = -5.53;
	c[4] = -13.95;
	c[5] = -22.75;
	c[6] = -36.15;
	c[7] = -53.24;
	c[8] = -45.27;

	d[1] = -8.91;
	d[2] = -3.95;
	d[3] = -8.77;
	d[4] = -13.90;
	d[5] = -23.74;
	d[6] = -28.83;
	d[7] = -14.03;
	d[8] = -7.98;




	category_number = -1;
	for (i=1; i<=category_total_number; i++)
	{
		if ( (skyclearness >= category_bounds[i]) && (skyclearness < category_bounds[i+1]) )
			category_number = i;
	}

	if (category_number == -1) {
		if (suppress_warnings==0)
		    fprintf(stderr, "ERROR: Model parameters out of range\n");
		print_error_sky();
		exit(1);
	}
		

	value = a[category_number] + b[category_number]*atm_preci_water  + c[category_number]*cos(sunzenith*M_PI/180) + 
	    d[category_number]*log(skybrightness);

	return(value);
}


/* direct normal efficacy model, according to PEREZ */

double direct_n_effi_PEREZ()

{
double 	value;
double 	category_bounds[10], a[10], b[10], c[10], d[10];
int   	category_total_number, category_number, i;


if ((skyclearness<skyclearinf || skyclearness>skyclearsup || skybrightness<skybriginf || skybrightness>skybrigsup) && suppress_warnings==0)
   fprintf(stderr, "Warning: skyclearness or skybrightness out of range in function direct_n_effi_PEREZ \n");


/* initialize category bounds (clearness index bounds) */

category_total_number = 8;

category_bounds[1] = 1;
category_bounds[2] = 1.065;
category_bounds[3] = 1.230;
category_bounds[4] = 1.500;
category_bounds[5] = 1.950;
category_bounds[6] = 2.800;
category_bounds[7] = 4.500;
category_bounds[8] = 6.200;
category_bounds[9] = 12.1;


/* initialize model coefficients */
a[1] = 57.20;
a[2] = 98.99;
a[3] = 109.83;
a[4] = 110.34;
a[5] = 106.36;
a[6] = 107.19;
a[7] = 105.75;
a[8] = 101.18;

b[1] = -4.55;
b[2] = -3.46;
b[3] = -4.90;
b[4] = -5.84;
b[5] = -3.97;
b[6] = -1.25;
b[7] = 0.77;
b[8] = 1.58;

c[1] = -2.98;
c[2] = -1.21;
c[3] = -1.71;
c[4] = -1.99;
c[5] = -1.75;
c[6] = -1.51;
c[7] = -1.26;
c[8] = -1.10;

d[1] = 117.12;
d[2] = 12.38;
d[3] = -8.81;
d[4] = -4.56;
d[5] = -6.16;
d[6] = -26.73;
d[7] = -34.44;
d[8] = -8.29;



for (i=1; i<=category_total_number; i++)
{
if ( (skyclearness >= category_bounds[i]) && (skyclearness < category_bounds[i+1]) )
category_number = i;
}

value = a[category_number] + b[category_number]*atm_preci_water  + c[category_number]*exp(5.73*sunzenith*M_PI/180 - 5) +  d[category_number]*skybrightness;

if (value < 0) value = 0;

return(value);
}


/*check the range of epsilon and delta indexes of the perez parametrization*/
void check_parametrization()
{
if (skyclearness<skyclearinf || skyclearness>skyclearsup || skybrightness<skybriginf || skybrightness>skybrigsup)
		{

/*  limit sky clearness or sky brightness, 2009 11 13 by J. Wienold */
		if (skyclearness<skyclearinf){
			skyclearness=skyclearinf;
			if (suppress_warnings==0)
			    fprintf(stderr,"Range warning: sky clearness too low (%lf)\n", skyclearness);
		}
		if (skyclearness>skyclearsup){
			skyclearness=skyclearsup-0.1;
			if (suppress_warnings==0)
			    fprintf(stderr,"Range warning: sky clearness too high (%lf)\n", skyclearness);
		}
		if (skybrightness<skybriginf){
			skybrightness=skybriginf;
			if (suppress_warnings==0)
			    fprintf(stderr,"Range warning: sky brightness too low (%lf)\n", skybrightness);
		}
		if (skybrightness>skybrigsup){
			skybrightness=skybrigsup;
			if (suppress_warnings==0)
			    fprintf(stderr,"Range warning: sky brightness too high (%lf)\n", skybrightness);
		}

	return;	}
	else return;
}


/* validity of the direct and diffuse components */
void 	check_illuminances()
{
	if (directilluminance < 0) {
		fprintf(stderr,"WARNING: direct illuminance < 0. Using 0.0\n");
		directilluminance = 0.0;
	}
	if (diffuseilluminance < 0) {
		fprintf(stderr,"WARNING: diffuse illuminance < 0. Using 0.0\n");
		diffuseilluminance = 0.0;
	}
	if (directilluminance > solar_constant_l*1000.0) {
		fprintf(stderr,"ERROR: direct illuminance exceeds solar constant\n");
		exit(1);
	}
}


void 	check_irradiances()
{
	if (directirradiance < 0) {
		fprintf(stderr,"WARNING: direct irradiance < 0. Using 0.0\n");
		directirradiance = 0.0;
	}
	if (diffuseirradiance < 0) {
		fprintf(stderr,"WARNING: diffuse irradiance < 0. Using 0.0\n");
		diffuseirradiance = 0.0;
	}
	if (directirradiance > solar_constant_e) {
		fprintf(stderr,"ERROR: direct irradiance exceeds solar constant\n");
		exit(1);
	}
}
	


/* Perez sky's brightness */
double sky_brightness()
{
double value;

value = diffuseirradiance * air_mass() / ( solar_constant_e*get_eccentricity());

return(value);
}


/* Perez sky's clearness */
double sky_clearness()
{
	double value;

	value = ( (diffuseirradiance + directirradiance)/(diffuseirradiance) + 1.041*sunzenith*M_PI/180*sunzenith*M_PI/180*sunzenith*M_PI/180 ) / (1 + 1.041*sunzenith*M_PI/180*sunzenith*M_PI/180*sunzenith*M_PI/180) ;

	return(value);
}



/* diffus horizontal irradiance from Perez sky's brightness */
double diffuse_irradiance_from_sky_brightness()
{
	double value;

	value = skybrightness / air_mass() * ( solar_constant_e*get_eccentricity());

	return(value);
}


/* direct normal irradiance from Perez sky's clearness */
double direct_irradiance_from_sky_clearness()
{
	double value;

	value = diffuse_irradiance_from_sky_brightness();
	value = value * ( (skyclearness-1) * (1+1.041*sunzenith*M_PI/180*sunzenith*M_PI/180*sunzenith*M_PI/180) );

	return(value);
}


void illu_to_irra_index()
{
double	test1=0.1, test2=0.1, d_eff;
int	counter=0;	

diffuseirradiance = diffuseilluminance*solar_constant_e/(solar_constant_l*1000);
directirradiance = directilluminance*solar_constant_e/(solar_constant_l*1000);
skyclearness =  sky_clearness(); 
skybrightness = sky_brightness();
check_parametrization();
	
while ( ((fabs(diffuseirradiance-test1)>10) || (fabs(directirradiance-test2)>10)
		|| skyclearness>skyclearinf || skyclearness<skyclearsup 
		|| skybrightness>skybriginf || skybrightness<skybrigsup )
		 && !(counter==5) )
	{

	test1=diffuseirradiance;
	test2=directirradiance;	
	counter++;
	
	diffuseirradiance = diffuseilluminance/glob_h_diffuse_effi_PEREZ();
	d_eff = direct_n_effi_PEREZ();
	if (d_eff < 0.1)
		directirradiance = 0;
	else
		directirradiance = directilluminance/d_eff;
	
	skybrightness = sky_brightness();
	skyclearness =  sky_clearness();
	check_parametrization();
	
	}


return;
}		

static int get_numlin(float epsilon)
{
	if (epsilon < 1.065)
		return 0;
	else if (epsilon < 1.230) 
		return 1;
	else if (epsilon < 1.500)
		return 2;
	else if (epsilon < 1.950)
		return 3;
	else if (epsilon < 2.800)
		return 4;
	else if (epsilon < 4.500)
		return 5;
	else if (epsilon < 6.200)
		return 6;
	return 7;
}

/* sky luminance perez model */
double calc_rel_lum_perez(double dzeta,double gamma,double Z,double epsilon,double Delta,float coeff_perez[])
{
	float x[5][4];
	int i,j,num_lin;
	double c_perez[5];

	if ( (epsilon <  skyclearinf) || (epsilon >= skyclearsup) )
	{
		fprintf(stderr,"Epsilon out of range in function calc_rel_lum_perez!\n");
		exit(1);
	}

	/* correction de modele de Perez solar energy ...*/
	if ( (epsilon > 1.065) && (epsilon < 2.8) )
	{
		if ( Delta < 0.2 ) Delta = 0.2;
	}

	num_lin = get_numlin(epsilon);

	for (i=0;i<5;i++)
		for (j=0;j<4;j++)
		{
			x[i][j] = *(coeff_perez + 20*num_lin + 4*i +j);
			/* printf("x %d %d vaut %f\n",i,j,x[i][j]); */
		}


	if (num_lin)
	{
		for (i=0;i<5;i++)
			c_perez[i] = x[i][0] + x[i][1]*Z + Delta * (x[i][2] + x[i][3]*Z);
	}
	else
	{
		c_perez[0] = x[0][0] + x[0][1]*Z + Delta * (x[0][2] + x[0][3]*Z);
		c_perez[1] = x[1][0] + x[1][1]*Z + Delta * (x[1][2] + x[1][3]*Z);
		c_perez[4] = x[4][0] + x[4][1]*Z + Delta * (x[4][2] + x[4][3]*Z);
		c_perez[2] = exp( pow(Delta*(x[2][0]+x[2][1]*Z),x[2][2])) - x[2][3];
		c_perez[3] = -exp( Delta*(x[3][0]+x[3][1]*Z) )+x[3][2]+Delta*x[3][3];
	}


	return (1 + c_perez[0]*exp(c_perez[1]/cos(dzeta)) ) *
	    (1 + c_perez[2]*exp(c_perez[3]*gamma) +
	    c_perez[4]*cos(gamma)*cos(gamma) );
}



/* coefficients for the sky luminance perez model */
void coeff_lum_perez(double Z, double epsilon, double Delta, float coeff_perez[])
{
	float x[5][4];
	int i,j,num_lin;

	if ( (epsilon <  skyclearinf) || (epsilon >= skyclearsup) )
	{
		fprintf(stderr,"Epsilon out of range in function calc_rel_lum_perez !\n");
		exit(1);
	}

	/* correction du modele de Perez solar energy ...*/
	if ( (epsilon > 1.065) && (epsilon < 2.8) )
	{
		if ( Delta < 0.2 ) Delta = 0.2;
	}

	num_lin = get_numlin(epsilon);

	//fprintf(stderr,"numlin %d\n", num_lin);

	for (i=0;i<5;i++)
		for (j=0;j<4;j++)
		{
			x[i][j] = *(coeff_perez + 20*num_lin + 4*i +j);
			/* printf("x %d %d vaut %f\n",i,j,x[i][j]); */
		}


	if (num_lin)
	{
		for (i=0;i<5;i++)
			*(c_perez+i) = x[i][0] + x[i][1]*Z + Delta * (x[i][2] + x[i][3]*Z);

	}
	else
	{
		*(c_perez+0) = x[0][0] + x[0][1]*Z + Delta * (x[0][2] + x[0][3]*Z);
		*(c_perez+1) = x[1][0] + x[1][1]*Z + Delta * (x[1][2] + x[1][3]*Z);
		*(c_perez+4) = x[4][0] + x[4][1]*Z + Delta * (x[4][2] + x[4][3]*Z);
		*(c_perez+2) = exp( pow(Delta*(x[2][0]+x[2][1]*Z),x[2][2])) - x[2][3];
		*(c_perez+3) = -exp( Delta*(x[3][0]+x[3][1]*Z) )+x[3][2]+Delta*x[3][3];


	}


	return;
}


/* degrees into radians */
double radians(double degres)
{
	return degres*M_PI/180.0;
}

/* radian into degrees */
double degres(double radians)
{
	return radians/M_PI*180.0;
}

/* calculation of the angles dzeta and gamma */
void theta_phi_to_dzeta_gamma(double theta,double phi,double *dzeta,double *gamma, double Z)
{
	*dzeta = theta; /* dzeta = phi */
	if ( (cos(Z)*cos(theta)+sin(Z)*sin(theta)*cos(phi)) > 1 && (cos(Z)*cos(theta)+sin(Z)*sin(theta)*cos(phi) < 1.1 ) )
		*gamma = 0;
	else if ( (cos(Z)*cos(theta)+sin(Z)*sin(theta)*cos(phi)) > 1.1 )
	{
		printf("error in calculation of gamma (angle between point and sun");
		exit(3);
	}
	else
		*gamma = acos(cos(Z)*cos(theta)+sin(Z)*sin(theta)*cos(phi));
}



/********************************************************************************/
/*	Fonction: integ_lv							*/
/*										*/
/*	In: float *lv,*theta							*/
/*	    int sun_pos								*/
/*										*/
/*	Out: double								*/
/*										*/
/*	Update: 29/08/93							*/
/*										*/
/*	Rem: 									*/
/*										*/
/*	But: calcul l'integrale de luminance relative sans la dir. du soleil	*/
/*										*/
/********************************************************************************/
double integ_lv(float *lv,float *theta)
{
	int i;
	double buffer=0.0;

	for (i=0;i<145;i++)
		buffer += (*(lv+i))*cos(radians(*(theta+i)));

	return buffer*2*M_PI/144;

}






/* enter day number(double), return E0 = square(R0/R):  eccentricity correction factor  */

double get_eccentricity()
{
	double day_angle;
	double E0;

	day_angle  = 2*M_PI*(daynumber - 1)/365;
	E0         = 1.00011+0.034221*cos(day_angle)+0.00128*sin(day_angle)+
	    0.000719*cos(2*day_angle)+0.000077*sin(2*day_angle);

	return (E0);

}


/* enter sunzenith angle (degrees) return relative air mass (double) */
double 	air_mass()
{
double	m;

if (sunzenith>90)
	{
	fprintf(stderr, "Solar zenith angle larger than 90 degrees in function air_mass()\n");
	exit(1);
	}
	
m = 1/( cos(sunzenith*M_PI/180)+0.15*exp( log(93.885-sunzenith)*(-1.253) ) );
return(m);
}



