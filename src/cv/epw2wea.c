#ifndef lint
static const char RCSid[] = "$Id: epw2wea.c,v 2.5 2024/08/02 21:18:40 greg Exp $";
#endif
/*  Copyright (c) 2003
 *  National Research Council Canada
 *  written by Christoph Reinhart
 */

/* epw2wea: daylight analysis subprogram of DAYSIM */
/* Program converts EnergyPlus weather format (*.ppw) into DAYSIM format (*.wea) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main( int argc, char  *argv[])
{
	FILE *EPW_FILE = NULL;
	FILE* WEA_FILE = stdout;
	int year, month,day, hour_in,minute=60,i;
	int minute_message=1;
	int get_ac = 0;
	float dir_norm_rad, dif_or_rad,dummy_float;
	float aod;
	int cc;
    char keyword[2000]="";
	char minute_string[2000]="";
	char *epw_file = NULL;
	char *wea_file = NULL;
	char city[200] ="";
	char country[200] ="";
	char latitude[200] ="",longitude[200] ="",time_zone[200] ="",elevation[200] ="";

	for ( ; argc > 2 && argv[1][0] == '-'; argc--, argv++)
		if (argv[1][1] == 'a') {
			get_ac = !get_ac;
		} else {
			fprintf(stderr, "epw2wea: unknown option: %s\n", argv[1]);
			exit(1);
		}

	if ((argc < 2) | (argc > 3))
	{
		fprintf(stderr,"epw2wea: FATAL ERROR - wrong number of arguments\n");
		fprintf(stderr,"Usage: epw2wea [-a] file-name.epw [file-name.wea]\n");
		exit(1);
	}
	epw_file = argv[1];
	EPW_FILE=fopen(epw_file, "r");
	if (!EPW_FILE) {
		fprintf(stderr, "%s: cannot open for reading\n", epw_file);
		exit(1);
	}
	if (argc == 3) {
		wea_file = argv[2];
		WEA_FILE=fopen(wea_file, "w");
		if (!WEA_FILE) {
			fprintf(stderr, "%s: cannot open for writing\n", wea_file);
			exit(1);
		}
	}
	fscanf(EPW_FILE,"%[^,]s",keyword);
	if( !strcmp(keyword,"LOCATION") ){
	fscanf(EPW_FILE,",%[^,]s",city);
	fscanf(EPW_FILE,",%[^,]s",country);
	fscanf(EPW_FILE,",%[^,]s",country);
	sprintf(keyword,"place %s_%s\n",city,country);
	printf("%s",keyword);
	fprintf(WEA_FILE,"%s",keyword);

	fscanf(EPW_FILE,",%[^,]s",country);
	fscanf(EPW_FILE,",%[^,]s",country);
	fscanf(EPW_FILE,",%[^,]s",latitude);
	printf("latitude %s\n",latitude);
	fprintf(WEA_FILE,"latitude %s\n",latitude);
	fscanf(EPW_FILE,",%[^,]s",longitude);

	printf("longitude %.2f\n",-1.0*atof(longitude));
	fprintf(WEA_FILE,"longitude %.2f\n",-1.0*atof(longitude));
	fscanf(EPW_FILE,",%[^,]s",time_zone);
	printf("time_zone %.2f\n",-15.0*atof(time_zone));
	fprintf(WEA_FILE,"time_zone %.0f\n",-15.0*atoi(time_zone));
	fscanf(EPW_FILE,",%s[^\n]",elevation);
	printf("site_elevation %s\nweather_data_file_units 1\n",elevation);
	fprintf(WEA_FILE,"site_elevation %s\nweather_data_file_units 1\n",elevation);

	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		
	/* read in time step interval */
	fscanf(EPW_FILE,"%[^,]s",keyword);
	fscanf(EPW_FILE,",%[^,]s",keyword);
	fscanf(EPW_FILE,",%[^,]s",minute_string);
	minute=atoi(minute_string);
	if(minute==1)   /* one measurement per hour equals a 60 minute time step */
		minute=60;
	fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
    
	while( EOF != fscanf(EPW_FILE,"%d,%d,%d,%d",&year,&month,&day, &hour_in)){
		
		fprintf(WEA_FILE,"%d %d %.3f ",month,day,hour_in*1.0-minute*(0.5/60));
		
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%[^,]s",city);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);
		fscanf(EPW_FILE,",%f",&dummy_float);

		fscanf(EPW_FILE,",%f,%f",&dir_norm_rad, &dif_or_rad);
		fprintf(WEA_FILE,"%.0f %.0f",dir_norm_rad, dif_or_rad);
		if (get_ac){
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%d",&cc);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&aod);

			fprintf(WEA_FILE," %.3f %.1f",aod, cc/10.);
			
		}

		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fprintf(WEA_FILE,"\n");
	}

}else{
		fprintf(stderr,"epw2wea: FATAL ERROR - this does not seem to be an epw file \n");exit(1);
}

fclose(EPW_FILE);
fclose(WEA_FILE);
return 0;
}


