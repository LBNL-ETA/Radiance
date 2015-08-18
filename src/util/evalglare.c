#ifndef lint
static const char RCSid[] = "$Id: evalglare.c,v 2.2 2015/08/18 15:02:53 greg Exp $";
#endif
/* EVALGLARE V1.17
 * Evalglare Software License, Version 1.0
 *
 * Copyright (c) 1995 - 2015 Fraunhofer ISE, EPFL.
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, WITHOUT
 * modification, are permitted provided that the following all conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgments:
 *             "This product includes the evalglare software,
                developed at Fraunhofer ISE by J. Wienold" and 
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Evalglare," and "Fraunhofer ISE" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact jan.wienold@epfl.ch
 *
 * 5. Products derived from this software may not be called "evalglare",
 *       nor may "evalglare" appear in their name, without prior written
 *       permission of Fraunhofer ISE.
 *
 * Redistribution and use in source and binary forms, WITH 
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *
 * conditions 1.-5. apply
 * 
 * 6. In order to ensure scientific correctness, any modification in source code imply fulfilling all following comditions:
 *     a) A written permission from Fraunhofer ISE. For written permission, please contact jan.wienold@ise.fraunhofer.de.
 *     b) The permission can be applied via email and must contain the applied changes in source code and at least two example calculations, 
 *        comparing the results of the modified version and the current version of evalglare. 
 *     b) Any redistribution of a modified version of evalglare must contain following note:
 *        "This software uses a modified version of the source code of evalglare."
 *
 *
 *    
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Fraunhofer ISE  OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/).
 *
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2009 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* evalglare.c, v0.2 2005/08/21 18:00:00 wienold */
/* evalglare.c, v0.3 2006/06/01 16:20:00 wienold
   changes to the v02 version:
   -fix problem with non-square pictures
   -set luminance values in the hemisphere behind the line of sight to 0
   -fix error in vcp calculation */
/* evalglare.c, v0.4 2006/11/13 15:10:00 wienold
   changes to the v03 version:
   -fix problem with tabulator in picture header "cannot read view"
   -fix missing pixels in task area
   evalglare.c, v0.5 2006/11/29 9:30 wienold
   changes to the v04 version:
   - fix problem with glare sources at the edge of the picture
*/
/* evalglare.c, v0.6 2006/11/29 9:30 wienold
   changes to the v05 version:
   - add luminance restriction parameters in order to take into account restrictions of measuring systems
*/

/* evalglare.c, v0.7 2007/07/11 18:00 wienold
   changes to the v06 version:
   - add external provision of vertical illuminance
*/

/* evalglare.c, v0.8 2007/12/04 1:35 wienold
   changes to the v07 version:
   - limit dgp to max 1.0
   - fill up cutted pictures with last known values
   - add second detailed output
*/

/* evalglare.c, v0.9 2008/07/02  wienold
   changes to the v08 version:
   - add version number in picheader 
   - add option for showing version numer
*/
/* evalglare.c, v0.9a 2008/09/20  wienold
   changes to the v09 version:
   - reduce fix value threshold from 500 to 100 
*/

/* evalglare.c, v0.9b 2008/11/12  wienold
   changes to the v09a version:
   - check view type vta, if wrong than exit 
*/
/* evalglare.c, v0.9c 2009/03/31  wienold
   changes to the v09b version:
   - peak extraction is default now (-y) , for deactivation use -x
*/

/* evalglare.c, v0.9d 2009/06/24  wienold
   changes to the v09c version:
   - fixed memory problem while using strcat
*/
/* evalglare.c, v0.9e 2009/10/10  wienold
   changes to the v09d version:
   - fixed problem while reading .pic of windows version
*/
/* evalglare.c, v0.9f 2009/10/21  wienold
   changes to the v09e version:
   - piping of input pictures possible, allow also vtv and vtc viewtyps
*/

/* evalglare.c, v0.9g 2009/10/21  wienold
   changes to the v09f version:
   - modified pictool.c
    - added -V (calc only vertical illuminance)
*/
/* evalglare.c, v0.9h 2011/10/10  wienold
   changes to the v09g version:
   - include disability glare for age factor of 1
   - 
*/
/* evalglare.c, v0.9i 2011/10/17  wienold
   changes to the v09h version:
   - M option: Correct wrong (measured) luminance images by feeding Ev and a luminance value, which should be replaced  
    
*/
/* evalglare.c, v1.0 2012/02/08  wienold
   changes to the v09h version:
   - include all view types
   - check view option in header
   - add view options in command line
   - option for cutting out GUTH field of view   
    
*/

/* evalglare.c, v1.02 2012/03/01  wienold,reetz,grobe
   changes to the v1.0 version:
   - fixed buffer overflow for string variable version[40]
   - replaced various strings to handle release by #define
   - removed most unused variables
   - initialized variables
   - removed nested functions
   - compiles now with -ansi -pedantic
  
*/

/* evalglare.c, v1.03 2012/04/17  wienold
   - include low light correction
*/


/* evalglare.c, v1.04 2012/04/23  wienold
   - remove bug for gen_dgp_profile output
*/

/* evalglare.c, v1.05 2012/05/29  wienold
   - remove variable overflow of low-light-correction for large Ev
*/
/* evalglare.c, v1.06 2012/05/29  wienold
   - initiate low-light-correction-variable
*/
/* evalglare.c, v1.07 2012/05/29  wienold
   - remove edge pixels from evaluation, when center of pixel is behind view (stability)
*/
 
/* evalglare.c, v1.08 2012/09/09  wienold
   - add direction vector in detailed output for each glare source
   - include age correction 
*/
/* evalglare.c, v1.09 2012/09/09  wienold
   - faster calculation for the use of gen_dgp_profile: no vertical illuminance calculation, only dgp is calculated, second scan + split is deactivated, when no pixel >15000 is found 
*/
/* evalglare.c, v1.10 2012/09/09  wienold
   - faster calculation for the use of gen_dgp_profile: no vertical illuminance calculation, only dgp is calculated, second scan + split is deactivated, when no pixel >15000 is found 
*/
/* evalglare.c, v1.11 2013/01/17  wienold
   - fix output bug of dgp, when using -i or -I
   */
   
/* evalglare.c, v1.12 2013/10/31  wienold
   - include CIE equation for disability glare, Stiles-Holladay
   */
/* evalglare.c, v1.13 2014/04/06  wienold
   - remove bug: initialize Lveil_cie_sum
   used for the CIE equation for disability glare
   */

/* evalglare.c, v1.14 buggy changes... removed...
   */
/* evalglare.c, v1.15 add option for uniform colored glare sources
   */
/* evalglare.c, v1.16 2015/05/05   remove bugs: background luminance is now calculated not cos-weighted any more, to switch on cos-weighting use -n option
                                   calculation of the background luminance is now based on the real solid angle and not hard coded 2*PI any more
				   fast calculation mode: total solid angle =2*PI (before PI)
   */
   
/* evalglare.c, v1.17 2015/07/15 add option for calculating band luminance -B angle_of_band
   remove of age factor due to non proven statistical evidence 
   */

#define EVALGLARE

#define PROGNAME "evalglare"
#define VERSION "1.17 release 15.07.2015 by EPFL, J.Wienold"
#define RELEASENAME PROGNAME " " VERSION

#include "rtio.h"
#include "platform.h"
#include "pictool.h"
#include <math.h>
#include <string.h>
char *progname;

/* subroutine to add a pixel to a glare source */
void add_pixel_to_gs(pict * p, int x, int y, int gsn)
{
	double old_av_posx, old_av_posy, old_av_lum, old_omega, act_omega,
		new_omega, act_lum;


	pict_get_npix(p, gsn) = pict_get_npix(p, gsn) + 1;
	old_av_posx = pict_get_av_posx(p, gsn);
	old_av_posy = pict_get_av_posy(p, gsn);
	old_av_lum = pict_get_av_lum(p, gsn);
	old_omega = pict_get_av_omega(p, gsn);

	act_omega = pict_get_omega(p, x, y);
	act_lum = luminance(pict_get_color(p, x, y));
	new_omega = old_omega + act_omega;
	pict_get_av_posx(p, gsn) =
		(old_av_posx * old_omega + x * act_omega) / new_omega;
	pict_get_av_posy(p, gsn) =
		(old_av_posy * old_omega + y * act_omega) / new_omega;
	if (isnan((pict_get_av_posx(p, gsn))))
		fprintf(stderr,"error in add_pixel_to_gs %d %d %f %f %f %f\n", x, y, old_av_posy, old_omega, act_omega, new_omega);

	pict_get_av_lum(p, gsn) =
		(old_av_lum * old_omega + act_lum * act_omega) / new_omega;
	pict_get_av_omega(p, gsn) = new_omega;
	pict_get_gsn(p, x, y) = gsn;
	if (act_lum < pict_get_lum_min(p, gsn)) {
		pict_get_lum_min(p, gsn) = act_lum;
	}
	if (act_lum > pict_get_lum_max(p, gsn)) {
		pict_get_lum_max(p, gsn) = act_lum;
	}

 /*    printf("gsn,x,y,av_posx,av_posy,av_lum  %i %f  %f %f %f %f\n",gsn,x,y,pict_get_av_posx(p, gsn),pict_get_av_posy(p, gsn),pict_get_av_lum(p, gsn)); */

}

/* subroutine for peak extraction */
int
find_split(pict * p, int x, int y, double r, int i_split_start,
		   int i_split_end)
{
	int i_find_split, x_min, x_max, y_min, y_max, ix, iy, iix, iiy, dx, dy,
		out_r;
	double r_actual;

	i_find_split = 0;

	x_min = 0;
	y_min = 0;
	x_max = pict_get_ysize(p) - 1;

	y_max = pict_get_ysize(p) - 1;

	for (iiy = 1; iiy <= 2; iiy++) {
		dy = iiy * 2 - 3;
		if (dy == -1) {
			iy = y;
		} else {
			iy = y + 1;
		}
		while (iy <= y_max && iy >= y_min) {
			out_r = 0;
			for (iix = 1; iix <= 2; iix++) {
				dx = iix * 2 - 3;
				if (dx == -1) {
					ix = x;
				} else {
					ix = x + 1;
				}
				while (ix <= x_max && ix >= x_min && iy >= y_min) {

					r_actual =
						acos(DOT(pict_get_cached_dir(p, x, y),
								 pict_get_cached_dir(p, ix, iy))) * 2;
					if (r_actual <= r) {
						out_r = 1;
						if (pict_get_gsn(p, ix, iy) >= i_split_start
							&& pict_get_gsn(p, ix, iy) <= i_split_end) {
							i_find_split = pict_get_gsn(p, ix, iy);
						}
					} else {
						ix = -99;
					}
					ix = ix + dx;
				}
			}
			if (out_r == 0) {
				iy = -99;
			}
			iy = iy + dy;

		}
	}



	return i_find_split;
}

/*
static	int	icomp(int* a,int* b)
{
	if (*a < *b)
		return -1;
	if (*a > *b)
		return 1;
	return 0;
}

*/

/* subroutine to find nearby glare source pixels */

int
find_near_pgs(pict * p, int x, int y, float r, int act_gsn, int max_gsn,
			  int min_gsn)
{
	int dx, dy, i_near_gs, xx, yy, x_min, x_max, y_min, y_max, ix, iy, iix,
		iiy, old_gsn, new_gsn, find_gsn, change, stop_y_search,
		stop_x_search;
	double r_actual;
	int ixm[3];

	i_near_gs = 0;
	stop_y_search = 0;
	stop_x_search = 0;
	x_min = 0;
	y_min = 0;
	if (act_gsn == 0) {
		x_max = x;
	} else {
		x_max = pict_get_xsize(p) - 1;
	}

	y_max = pict_get_ysize(p) - 1;

	old_gsn = pict_get_gsn(p, x, y);
	new_gsn = old_gsn;
	change = 0;
	if (act_gsn > 0) {
		i_near_gs = pict_get_gsn(p, x, y);
	}
	for (iiy = 1; iiy <= 2; iiy++) {
		dy = iiy * 2 - 3;
		if (dy == -1) {
			iy = y;
		} else {
			iy = y + 1;
		}
		ixm[1] = x;
		ixm[2] = x;
		stop_y_search = 0;


		while (iy <= y_max && iy >= y_min) {
			for (iix = 1; iix <= 2; iix++) {
				dx = iix * 2 - 3;
				ix = (ixm[1] + ixm[2]) / 2;
				stop_x_search = 0;
				while (ix <= x_max && ix >= x_min && stop_x_search == 0
					   && stop_y_search == 0) {
/*        printf(" dx,act_gsn, x,y,x_max, x_min, ix ,iy , ymax,ymin %i %i  %i %i %i  %i %i %i %i %i\n",dx,act_gsn,x,y,x_max,x_min,ix,iy,y_max,y_min);*/
					r_actual =
						acos(DOT(pict_get_cached_dir(p, x, y),
								 pict_get_cached_dir(p, ix, iy))) * 2;
					if (r_actual <= r) {
						if (pict_get_gsn(p, ix, iy) > 0) {
							if (act_gsn == 0) {
								i_near_gs = pict_get_gsn(p, ix, iy);
								stop_x_search = 1;
								stop_y_search = 1;
							} else {
								find_gsn = pict_get_gsn(p, ix, iy);
								if (pict_get_av_omega(p, old_gsn) <
									pict_get_av_omega(p, find_gsn)
									&& pict_get_av_omega(p,
														 find_gsn) >
									pict_get_av_omega(p, new_gsn)
									&& find_gsn >= min_gsn
									&& find_gsn <= max_gsn) {
									/*  other primary glare source found with larger solid angle -> add together all */
									new_gsn = find_gsn;
									change = 1;
									stop_x_search = 1;
									stop_y_search = 1;
								}
							}
						}
					} else {
						ixm[iix] = ix - dx;
						stop_x_search = 1;
					}
					ix = ix + dx;
				}
			}
			iy = iy + dy;
		}
	}
	if (change > 0) {
		pict_get_av_lum(p, old_gsn) = 0.;
		pict_get_av_omega(p, old_gsn) = 0.;
		pict_get_npix(p, old_gsn) = 0.;
		pict_get_lum_max(p, old_gsn) = 0.;


		i_near_gs = new_gsn;
/*						    printf(" changing gs no %i",old_gsn);
						    printf(" to %i\n",new_gsn);
	  */
		for (xx = 0; xx < pict_get_xsize(p); xx++)
			for (yy = 0; yy < pict_get_ysize(p); yy++) {
				if (pict_is_validpixel(p, x, y)) {
					if (pict_get_gsn(p, xx, yy) == old_gsn) {
						add_pixel_to_gs(p, xx, yy, new_gsn);
					}
				}
			}
	}

	return i_near_gs;


}

/* subroutine for calculation of task luminance */
double get_task_lum(pict * p, int x, int y, float r, int task_color)
{
	int x_min, x_max, y_min, y_max, ix, iy;
	double r_actual, av_lum, omega_sum, act_lum;


	x_max = pict_get_xsize(p) - 1;
	y_max = pict_get_ysize(p) - 1;
	x_min = 0;
	y_min = 0;



	av_lum = 0.0;
	omega_sum = 0.0;

	for (iy = y_min; iy <= y_max; iy++) {
		for (ix = x_min; ix <= x_max; ix++) {

/*			if (DOT(pict_get_cached_dir(p,ix,iy),p->view.vdir) < 0.0) 
				continue;*/
			r_actual =
				acos(DOT(pict_get_cached_dir(p, x, y),
						 pict_get_cached_dir(p, ix, iy))) * 2;
			act_lum = luminance(pict_get_color(p, ix, iy));

			if (r_actual <= r) {
				act_lum = luminance(pict_get_color(p, ix, iy));
				av_lum += pict_get_omega(p, ix, iy) * act_lum;
				omega_sum += pict_get_omega(p, ix, iy);
				if (task_color == 1) {
					pict_get_color(p, ix, iy)[RED] = 0.0;
					pict_get_color(p, ix, iy)[GRN] = 0.0;
					pict_get_color(p, ix, iy)[BLU] =
						act_lum / WHTEFFICACY / CIE_bf;
				}
			}
		}
	}

	av_lum = av_lum / omega_sum;
/*    printf("average luminance of task %f \n",av_lum);
      printf("solid angle of task %f \n",omega_sum);*/
	return av_lum;
}

/* subroutine for calculation of band luminance */
double get_band_lum(pict * p, float r, int task_color)
{
	int x_min, x_max, y_min, y_max, ix, iy, y_mid;
	double r_actual, av_lum, omega_sum, act_lum;


	x_max = pict_get_xsize(p) - 1;
	y_max = pict_get_ysize(p) - 1;
	x_min = 0;
	y_min = 0;
	y_mid = rint(y_max/2);



	av_lum = 0.0;
	omega_sum = 0.0;

	for (iy = y_min; iy <= y_max; iy++) {
		for (ix = x_min; ix <= x_max; ix++) {

/*			if (DOT(pict_get_cached_dir(p,ix,iy),p->view.vdir) < 0.0) 
				continue;*/
			r_actual =
				acos(DOT(pict_get_cached_dir(p, ix, y_mid),
						 pict_get_cached_dir(p, ix, iy))) ;
			act_lum = luminance(pict_get_color(p, ix, iy));

			if ((r_actual <= r) || (iy == y_mid) ) {
				act_lum = luminance(pict_get_color(p, ix, iy));
				av_lum += pict_get_omega(p, ix, iy) * act_lum;
				omega_sum += pict_get_omega(p, ix, iy);
				if (task_color == 1) {
					pict_get_color(p, ix, iy)[RED] = 0.0;
					pict_get_color(p, ix, iy)[GRN] = act_lum / WHTEFFICACY / CIE_gf;
					pict_get_color(p, ix, iy)[BLU] = 0.0;
				}
			}
		}
	}

	av_lum = av_lum / omega_sum;
/*    printf("average luminance of band %f \n",av_lum);*/
/*      printf("solid angle of band %f \n",omega_sum);*/
	return av_lum;
}





/* subroutine for coloring the glare sources 
     red is used only for explicit call of the subroutine with acol=0  , e.g. for disability glare
     -1: set to grey*/
int setglcolor(pict * p, int x, int y, int acol, int uniform_gs, double u_r, double u_g ,double u_b)
{
	int icol;
	double pcol[3][1000], act_lum, l;

        l=u_r+u_g+u_b ;
	
	pcol[RED][0] = 1.0 / CIE_rf;
	pcol[GRN][0] = 0.;
	pcol[BLU][0] = 0.;

	pcol[RED][1] = 0.;
	pcol[GRN][1] = 1.0 / CIE_gf;
	pcol[BLU][1] = 0.;

	pcol[RED][2] = 0.;
	pcol[GRN][2] = 0.;
	pcol[BLU][2] = 1 / CIE_bf;

	pcol[RED][3] = 1.0 / (1.0 - CIE_bf);
	pcol[GRN][3] = 1.0 / (1.0 - CIE_bf);
	pcol[BLU][3] = 0.;

	pcol[RED][4] = 1.0 / (1.0 - CIE_gf);
	pcol[GRN][4] = 0.;
	pcol[BLU][4] = 1.0 / (1.0 - CIE_gf);

	pcol[RED][5] = 0.;
	pcol[GRN][5] = 1.0 / (1.0 - CIE_rf);
	pcol[BLU][5] = 1.0 / (1.0 - CIE_rf);

	pcol[RED][6] = 0.;
	pcol[GRN][6] = 1.0 / (1.0 - CIE_rf);
	pcol[BLU][6] = 1.0 / (1.0 - CIE_rf);


	pcol[RED][999] = 1.0 / WHTEFFICACY;
	pcol[GRN][999] = 1.0 / WHTEFFICACY;
	pcol[BLU][999] = 1.0 / WHTEFFICACY;
 
 
 	pcol[RED][998] = u_r /(l* CIE_rf) ;
	pcol[GRN][998] = u_g /(l* CIE_gf);
	pcol[BLU][998] = u_b /(l* CIE_bf);
/*        printf("CIE_rf,CIE_gf,CIE_bf,l=%f %f %f %f\n",CIE_rf,CIE_gf,CIE_bf,l);*/
        icol = acol ;
	if ( acol == -1) {icol=999;
	                          }else{if (acol>0){icol = acol % 5 +1;
	                         }};
        if ( uniform_gs == 1) {icol=998;
	                          } ;
	
	act_lum = luminance(pict_get_color(p, x, y));

	pict_get_color(p, x, y)[RED] = pcol[RED][icol] * act_lum / WHTEFFICACY;
	pict_get_color(p, x, y)[GRN] = pcol[GRN][icol] * act_lum / WHTEFFICACY;
	pict_get_color(p, x, y)[BLU] = pcol[BLU][icol] * act_lum / WHTEFFICACY;

	return icol;
}

/* this is the smooting subroutine */
void add_secondary_gs(pict * p, int x, int y, float r, int act_gs, int uniform_gs, double u_r, double u_g ,double u_b)
{
	int x_min, x_max, y_min, y_max, ix, iy, icol;
	double r_actual, omega_gs, omega_total, om;



	omega_gs = 0.0;
	omega_total = 0.0;
	x_min = x - r;
	if (x_min < 0) {
		x_min = 0;
	}
	x_max = x + r;
	if (x_max > pict_get_xsize(p) - 1) {
		x_max = pict_get_xsize(p) - 2;
	}
	y_min = y - r;
	if (y_min < 0) {
		y_min = 0;
	}
	y_max = y + r;
	if (y_max > pict_get_ysize(p) - 1) {
		y_max = pict_get_ysize(p) - 2;
	}

	for (ix = x_min; ix <= x_max; ix++)
		for (iy = y_min; iy <= y_max; iy++) {
			r_actual = sqrt((x - ix) * (x - ix) + (y - iy) * (y - iy));
			if (r_actual <= r) {
				om = pict_get_omega(p, ix, iy);
				omega_total += om;
				if (pict_get_gsn(p, ix, iy) == act_gs
					&& pict_get_pgs(p, ix, iy) == 1) {
					omega_gs = omega_gs + 1 * om;
				}
			}
		}


	if (omega_gs / omega_total > 0.2) {
/* add pixel to gs */

		add_pixel_to_gs(p, x, y, act_gs);

/* color pixel of gs */

		icol = setglcolor(p, x, y, act_gs, uniform_gs, u_r, u_g , u_b);

	}
}

/* subroutine for removing a pixel from a glare source */
void split_pixel_from_gs(pict * p, int x, int y, int new_gsn, int uniform_gs, double u_r, double u_g , double u_b)
{
	int old_gsn, icol;
	double old_av_posx, old_av_posy, old_av_lum, old_omega, act_omega,
		new_omega, act_lum;


/* change existing gs */
	old_gsn = pict_get_gsn(p, x, y);
	pict_get_npix(p, old_gsn) = pict_get_npix(p, old_gsn) - 1;

	act_omega = pict_get_omega(p, x, y);
	old_av_posx = pict_get_av_posx(p, old_gsn);
	old_av_posy = pict_get_av_posy(p, old_gsn);
	old_omega = pict_get_av_omega(p, old_gsn);

	new_omega = old_omega - act_omega;
	pict_get_av_omega(p, old_gsn) = new_omega;
	pict_get_av_posx(p, old_gsn) =
		(old_av_posx * old_omega - x * act_omega) / new_omega;
	pict_get_av_posy(p, old_gsn) =
		(old_av_posy * old_omega - y * act_omega) / new_omega;


	act_lum = luminance(pict_get_color(p, x, y));
	old_av_lum = pict_get_av_lum(p, old_gsn);
	pict_get_av_lum(p, old_gsn) =
		(old_av_lum * old_omega - act_lum * act_omega) / new_omega;
	/* add pixel to new  gs */

	add_pixel_to_gs(p, x, y, new_gsn);

/* color pixel of new gs */

	icol = setglcolor(p, x, y, new_gsn, uniform_gs, u_r, u_g , u_b);


}

/* subroutine for the calculation of the position index */
float get_posindex(pict * p, float x, float y, int postype)
{
	float posindex;
	double teta, phi, sigma, tau, deg, d, s, r, fact;


	pict_get_vangle(p, x, y, p->view.vdir, p->view.vup, &phi);
	pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &teta);
	pict_get_sigma(p, x, y, p->view.vdir, p->view.vup, &sigma);
	pict_get_tau(p, x, y, p->view.vdir, p->view.vup, &tau);

/* return (phi+teta+2*3.1415927); */

	deg = 180 / 3.1415927;
	fact = 0.8;
	if (phi == 0) {
		phi = 0.00001;
	}
	if (sigma <= 0) {
		sigma = -sigma;
	}
	if (teta == 0) {
		teta = 0.0001;
	}
	tau = tau * deg;
	sigma = sigma * deg;

	posindex =
		exp((35.2 - 0.31889 * tau -
			 1.22 * exp(-2 * tau / 9)) / 1000 * sigma + (21 +
														 0.26667 * tau -
														 0.002963 * tau *
														 tau) / 100000 *
			sigma * sigma);
/* below line of sight, using Iwata model */
	if (phi < 0) {
		d = 1 / tan(phi);
		s = tan(teta) / tan(phi);
		r = sqrt(1 / d * 1 / d + s * s / d / d);
		if (r > 0.6)
			fact = 1.2;
		if (r > 3) {
			fact = 1.2;
			r = 3;
		}

		posindex = 1 + fact * r;
	}
	if (posindex > 16)
		posindex = 16;

	return posindex;

}



double get_upper_cut_2eyes(float teta)
{
double phi;

phi=pow(7.7458218+0.00057407915*teta-0.00021746318*teta*teta+8.5572726e-6*teta*teta*teta,2);

return phi;
}
double get_lower_cut_2eyes(float teta)
{
double phi;

phi=1/(-0.014699242-1.5541106e-5*teta+4.6898068e-6*teta*teta-5.1539687e-8*teta*teta*teta);

return phi;
}

double get_lower_cut_central(float teta)
{
double phi;

phi=(68.227109-2.9524084*teta+0.046674262*teta*teta)/(1-0.042317294*teta+0.00075698419*teta*teta-6.5364285e-7*teta*teta*teta);

if (teta>73) {
 phi=60;
 }

return phi;
}

/* subroutine for the cutting the total field of view */
void cut_view_1(pict*p)
{
int x,y;
double border,ang,teta,phi,phi2;
		  for(x=0;x<pict_get_xsize(p)-1;x++)
		   for(y=0;y<pict_get_ysize(p)-1;y++) {
			if (pict_get_hangle(p,x,y,p->view.vdir,p->view.vup,&ang)) {
				if (pict_is_validpixel(p, x, y)) {
                       			pict_get_vangle(p,x,y,p->view.vdir,p->view.vup,&phi2);
					pict_get_sigma(p,x,y,p->view.vdir,p->view.vup,&phi);
					pict_get_tau(p,x,y,p->view.vdir,p->view.vup,&teta);
                                        phi=phi*180/3.1415927;
                                        phi2=phi2*180/3.1415927;
                                        teta=teta*180/3.1415927;
					if (teta<0) {
					teta=-teta;
					          }
                                        if(phi2>0){
					         border=get_upper_cut_2eyes(teta);

						 if (phi>border) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
					          }else{ 
					         border=get_lower_cut_2eyes(180-teta); 
						 if (-phi<border && teta>135) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
						  }
					
					
					
			     }else{
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;

			     
			     }

/*						 printf("teta,phi,border=%f %f %f\n",teta,phi,border);*/
			     }
			     }


return;

}
/* subroutine for the cutting the field of view seen by both eyes*/
void cut_view_2(pict*p)
{

int x,y;
double border,ang,teta,phi,phi2;
		  for(x=0;x<pict_get_xsize(p)-1;x++)
		   for(y=0;y<pict_get_ysize(p)-1;y++) {
			if (pict_get_hangle(p,x,y,p->view.vdir,p->view.vup,&ang)) {
				if (pict_is_validpixel(p, x, y)) {
                       			pict_get_vangle(p,x,y,p->view.vdir,p->view.vup,&phi2);
					pict_get_hangle(p,x,y,p->view.vdir,p->view.vup,&teta);
					pict_get_sigma(p,x,y,p->view.vdir,p->view.vup,&phi);
					pict_get_tau(p,x,y,p->view.vdir,p->view.vup,&teta);
                                        phi=phi*180/3.1415927;
                                        phi2=phi2*180/3.1415927;
                                        teta=teta*180/3.1415927;
/*						 printf("x,y,teta,phi,border=%i %i %f %f %f\n",x,y,teta,phi,border);*/
					if (teta<0) {
					teta=-teta;
					          }
                                        if(phi2>0){
					         border=60;

						 if (phi>border) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
					          }else{ 
					         border=get_lower_cut_central(180-teta);
						 if (phi>border) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
						  }
					
					
					
			     }else{
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;

			     
			     }
			     }
			     }


return;

}

/* subroutine for the cutting the field of view seen by both eyes
   luminance is modified by position index - just for experiments - undocumented */
void cut_view_3(pict*p)
{

int x,y;
double border,ang,teta,phi,phi2,lum,newlum;
		  for(x=0;x<pict_get_xsize(p)-1;x++)
		   for(y=0;y<pict_get_ysize(p)-1;y++) {
			if (pict_get_hangle(p,x,y,p->view.vdir,p->view.vup,&ang)) {
			     if (DOT(pict_get_cached_dir(p,x,y),p->view.vdir) >= 0.0) {
                       			pict_get_vangle(p,x,y,p->view.vdir,p->view.vup,&phi2);
					pict_get_hangle(p,x,y,p->view.vdir,p->view.vup,&teta);
					pict_get_sigma(p,x,y,p->view.vdir,p->view.vup,&phi);
					pict_get_tau(p,x,y,p->view.vdir,p->view.vup,&teta);
                                        phi=phi*180/3.1415927;
                                        phi2=phi2*180/3.1415927;
                                        teta=teta*180/3.1415927;
                                        lum=luminance(pict_get_color(p,x,y));
					
					newlum=lum/get_posindex(p,x,y,0);

		                        		pict_get_color(p,x,y)[RED]=newlum/WHTEFFICACY;
		                        		pict_get_color(p,x,y)[GRN]=newlum/WHTEFFICACY;
		                        		pict_get_color(p,x,y)[BLU]=newlum/WHTEFFICACY;

					
					
					
/*						 printf("x,y,teta,phi,border=%i %i %f %f %f\n",x,y,teta,phi,border);*/
					if (teta<0) {
					teta=-teta;
					          }
                                        if(phi2>0){
					         border=60;

						 if (phi>border) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
					          }else{ 
					         border=get_lower_cut_central(180-teta);
						 if (phi>border) {
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;
						 }
						  }
					
					
					
			     }else{
		                        		pict_get_color(p,x,y)[RED]=0;
		                        		pict_get_color(p,x,y)[GRN]=0;
		                        		pict_get_color(p,x,y)[BLU]=0;

			     
			     }
			     }
			     }


return;

}

/* subroutine for the calculation of the daylight glare index */


float get_dgi(pict * p, float lum_backg, int igs, int posindex_2)
{
	float dgi, sum_glare, omega_s;
	int i;


	sum_glare = 0;
	omega_s = 0;

	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			omega_s =
				pict_get_av_omega(p, i) / get_posindex(p,pict_get_av_posx(p,i), pict_get_av_posy(p,i),posindex_2) /
				get_posindex(p, pict_get_av_posx(p, i), pict_get_av_posy(p, i), posindex_2);
			sum_glare += 0.478 * pow(pict_get_av_lum(p, i), 1.6) * pow(omega_s,0.8) / (lum_backg + 0.07 * pow(pict_get_av_omega(p, i),0.5) * pict_get_av_lum(p, i));
                    /*    printf("i,sum_glare %i %f\n",i,sum_glare); */
		}
	}
	dgi = 10 * log10(sum_glare);

	return dgi;

}

/* subroutine for the calculation of the daylight glare probability */
double
get_dgp(pict * p, double E_v, int igs, double a1, double a2, double a3,
		double a4, double a5, double c1, double c2, double c3,
		int posindex_2)
{
	double dgp;
	double sum_glare;
	int i;

	sum_glare = 0;
	if (igs > 0) {
		for (i = 0; i <= igs; i++) {

			if (pict_get_npix(p, i) > 0) {
				sum_glare +=
					pow(pict_get_av_lum(p, i),
						a1) / pow(get_posindex(p, pict_get_av_posx(p, i),
											   pict_get_av_posy(p, i),
											   posindex_2),
								  a4) * pow(pict_get_av_omega(p, i), a2);
/*			printf("i,sum_glare %i %f\n",i,sum_glare);*/					  
			}
		}
		dgp =
			c1 * pow(E_v,
					 a5) + c3 + c2 * log10(1 + sum_glare / pow(E_v, a3));
	} else {
		dgp = c3 + c1 * pow(E_v, a5);
	}

	if (dgp > 1) {
		dgp = 1;
	}
/*			printf("dgp  %f\n",dgp);*/					  
	return dgp;

}

/* subroutine for the calculation of the visual comfort probability */
float get_vcp(pict * p, double lum_a, int igs, int posindex_2)
{
	float vcp;
	double sum_glare, dgr;
	int i, i_glare;


	sum_glare = 0;
	i_glare = 0;
	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			i_glare = i_glare + 1;
			sum_glare +=
				(0.5 * pict_get_av_lum(p, i) *
				 (20.4 * pict_get_av_omega(p, i) +
				  1.52 * pow(pict_get_av_omega(p, i),
							 0.2) - 0.075)) / (get_posindex(p,
															pict_get_av_posx
															(p, i),
															pict_get_av_posy
															(p, i),
															posindex_2) *
											   pow(lum_a, 0.44));

		}
	}
	dgr = pow(sum_glare, pow(i_glare, -0.0914));

	vcp = 50 * erf((6.347 - 1.3227 * log(dgr)) / 1.414213562373) + 50;
	if (dgr > 750) {
		vcp = 0;
	}
	if (dgr < 20) {
		vcp = 100;
	}


	return vcp;

}

/* subroutine for the calculation of the unified glare rating */
float get_ugr(pict * p, double lum_backg, int igs, int posindex_2)
{
	float ugr;
	double sum_glare;
	int i, i_glare;


	sum_glare = 0;
	i_glare = 0;
	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			i_glare = i_glare + 1;
			sum_glare +=
				pow(pict_get_av_lum(p, i) /
					get_posindex(p, pict_get_av_posx(p, i),
								 pict_get_av_posy(p, i), posindex_2),
					2) * pict_get_av_omega(p, i);

		}
	}
	ugr = 8 * log10(0.25 / lum_backg * sum_glare);

	return ugr;

}


/* subroutine for the calculation of the disability glare according to Poynter */
float get_disability(pict * p, double lum_backg, int igs)
{
	float disab;
	double sum_glare, sigma, deg;
	int i, i_glare;


	sum_glare = 0;
	i_glare = 0;
	deg = 180 / 3.1415927;
	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			i_glare = i_glare + 1;
			pict_get_sigma(p, pict_get_av_posx(p, i),
						   pict_get_av_posy(p, i), p->view.vdir,
						   p->view.vup, &sigma);

			sum_glare +=
				pict_get_av_lum(p,
								i) * cos(sigma +
										 0.00000000001) *
				pict_get_av_omega(p, i)
				/ (deg * sigma + 0.00000000001);
			if (isnan(sum_glare)) {
				printf("sigma for %f %f\n", pict_get_av_posx(p, i), pict_get_av_posy(p, i));
				printf("omega for %f %f\n", pict_get_av_posx(p, i), pict_get_av_posy(p, i));
                printf("avlum for %f %f\n", pict_get_av_posx(p, i), pict_get_av_posy(p, i));
                printf("avlum for %f %f %f\n", pict_get_av_posx(p, i), pict_get_av_posy(p, i), sigma);
			}

		}
	}
	disab = 5 * sum_glare;

	return disab;

}






/* subroutine for the calculation of the cie glare index */

float
get_cgi(pict * p, double E_v, double E_v_dir, int igs, int posindex_2)
{
	float cgi;
	double sum_glare;
	int i, i_glare;


	sum_glare = 0;
	i_glare = 0;
	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			i_glare = i_glare + 1;
			sum_glare +=
				pow(pict_get_av_lum(p, i) /
					get_posindex(p, pict_get_av_posx(p, i),
								 pict_get_av_posy(p, i), posindex_2),
					2) * pict_get_av_omega(p, i);

		}
	}
	cgi = 8 * log10((2 * (1 + E_v_dir / 500) / E_v) * sum_glare);

	return cgi;

}


#ifdef	EVALGLARE


/* main program 
------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char **argv)
{
	#define CLINEMAX 4095 /* memory allocated for command line string */
	pict *p = pict_create();
	int     skip_second_scan,calcfast,age_corr,cut_view,cut_view_type,calc_vill, output, detail_out2, y1, fill, yfillmax, yfillmin,
		ext_vill, set_lum_max, set_lum_max2, task_color, i_splitstart,
		i_split, posindex_2, task_lum, checkfile, rval, i, i_max, x, y,x2,y2,
		igs, actual_igs, icol, xt, yt, change, before_igs, sgs, splithigh,uniform_gs,
		detail_out, posindex_picture, non_cos_lb, rx, ry, rmx, rmy,apply_disability,band_calc,band_color;
	double  lum_total_max,age_corr_factor,age,dgp_ext,dgp,low_light_corr,omega_cos_contr, setvalue, lum_ideal, E_v_contr, sigma,
		E_vl_ext, lum_max, new_lum_max, r_center,
		search_pix, a1, a2, a3, a4, a5, c3, c1, c2, r_split, max_angle,
		omegat, sang, E_v, E_v2, E_v_dir, avlum, act_lum, ang,
		l_max, lum_backg, lum_backg_cos, omega_sources, lum_sources,
		lum, lum_source,teta,Lveil_cie,Lveil_cie_sum,disability_thresh,u_r,u_g,u_b,band_angle,band_avlum;
	float lum_task, lum_thres, dgi,  vcp, cgi, ugr, limit, 
		abs_max, Lveil;
	char file_out[500], file_out2[500], version[500];
	char *cline;
	VIEW userview = STDVIEW;
	int gotuserview = 0;

	/*set required user view parameters to invalid values*/
        uniform_gs = 0;
	apply_disability = 0;
	disability_thresh = 0;
	Lveil_cie_sum=0.0;
	skip_second_scan=0;
	lum_total_max=0.0;
	calcfast=0;
	age_corr_factor = 1.0;
	age_corr = 0;
	age = 20.0;
	userview.horiz = 0;
	userview.vert = 0;
	userview.type = 0;
        dgp_ext = 0;
	E_vl_ext = 0.0;
	new_lum_max = 0.0;
	lum_max = 0.0;
	omegat = 0.0;
	yt = 0;
	xt = 0;
	yfillmin = 0;
	yfillmax = 0;
	cut_view = 0;
	cut_view_type = 0;
	setvalue = 2e09;
	omega_cos_contr = 0.0;
	lum_ideal = 0.0;
	max_angle = 0.2;
	lum_thres = 5.0;
	task_lum = 0;
	sgs = 0;
	splithigh = 1;
	detail_out = 0;
	detail_out2 = 0;
	posindex_picture = 0;
	checkfile = 0;
	ext_vill = 0;
	fill = 0;
	a1 = 2.0;
	a2 = 1.0;
	a3 = 1.87;
	a4 = 2.0;
	a5 = 1.0;
	c1 = 5.87e-05;
	c2 = 0.092;
	c3 = 0.159;
	non_cos_lb = 1;
	posindex_2 = 0;
	task_color = 0;
	limit = 50000.0;
	set_lum_max = 0;
	set_lum_max2 = 0;
	abs_max = 0;
	progname = argv[0];
	E_v_contr = 0.0;
	strcpy(version, "1.15 release 05.02.2015 by J.Wienold");
        low_light_corr=1.0;
 	output = 0;
	calc_vill = 0;
	band_avlum = -99;
	band_calc = 0;
/* command line for output picture*/

	cline = (char *) malloc(CLINEMAX+1);
	cline[0] = '\0';
	for (i = 0; i < argc; i++) {
/*       fprintf(stderr, "%d %d \n",i,strlen(argv[i]));*/
		if (strlen(cline)+strlen(argv[i])+strlen(RELEASENAME)+2 >=CLINEMAX) {
			exit (-1);
		}
		else {
			strcat(cline, argv[i]);
			strcat(cline, " ");
		}
	}
	strcat(cline, "\n");
	strcat(cline, RELEASENAME);
	strcat(cline, "\n");


	/* program options */

	for (i = 1; i < argc && argv[i][0] == '-'; i++) {
		/* expand arguments */
		while ((rval = expandarg(&argc, &argv, i)) > 0);
		if (rval < 0) {
			fprintf(stderr, "%s: cannot expand '%s'", argv[0], argv[i]);
			exit(1);
		}
		rval = getviewopt(&userview, argc - i, argv + i);
		if (rval >= 0) {
			i += rval;
			gotuserview++;
			continue;
		}
		switch (argv[i][1]) {
		case 'a':
			age = atof(argv[++i]);
			age_corr = 1;
                        printf("age factor not supported any more \n");
			exit(1);
                        break;
		case 'b':
			lum_thres = atof(argv[++i]);
			break;
		case 'c':
			checkfile = 1;
			strcpy(file_out, argv[++i]);
			break;
		case 'u':
			uniform_gs = 1;
			u_r = atof(argv[++i]);
			u_g = atof(argv[++i]);
			u_b = atof(argv[++i]);
			break;
		case 'r':
			max_angle = atof(argv[++i]);
			break;
		case 's':
			sgs = 1;
			break;
		case 'k':
			apply_disability = 1;
			disability_thresh = atof(argv[++i]);
			break;
			
		case 'p':
			posindex_picture = 1;
			break;

		case 'y':
			splithigh = 1;
			break;
		case 'x':
			splithigh = 0;
			break;
		case 'Y':
			splithigh = 1;
			limit = atof(argv[++i]);
			break;

		case 'i':
			ext_vill = 1;
			E_vl_ext = atof(argv[++i]);
			break;
		case 'I':
			ext_vill = 1;
			fill = 1;
			E_vl_ext = atof(argv[++i]);
			yfillmax = atoi(argv[++i]);
			yfillmin = atoi(argv[++i]);
			break;
		case 'd':
			detail_out = 1;
			break;
		case 'D':
			detail_out2 = 1;
			break;
		case 'm':
			set_lum_max = 1;
			lum_max = atof(argv[++i]);
			new_lum_max = atof(argv[++i]);
			strcpy(file_out2, argv[++i]);
/*			printf("max lum set to %f\n",new_lum_max);*/
			break;

		case 'M':
			set_lum_max2 = 1;
			lum_max = atof(argv[++i]);
			E_vl_ext = atof(argv[++i]);
			strcpy(file_out2, argv[++i]);
/*			printf("max lum set to %f\n",new_lum_max);*/
			break;
		case 'n':
			non_cos_lb = 0;
			break;

		case 't':
			task_lum = 1;
			xt = atoi(argv[++i]);
			yt = atoi(argv[++i]);
			omegat = atof(argv[++i]);
			task_color = 0;
			break;
		case 'T':
			task_lum = 1;
			xt = atoi(argv[++i]);
			yt = atoi(argv[++i]);
/*			omegat= DEG2RAD(atof(argv[++i]));*/
			omegat = atof(argv[++i]);
			task_color = 1;
			break;
		case 'B':
			band_calc = 1;
			band_angle = atof(argv[++i]);
			band_color = 1;
			break;


		case 'w':
			a1 = atof(argv[++i]);
			a2 = atof(argv[++i]);
			a3 = atof(argv[++i]);
			a4 = atof(argv[++i]);
			a5 = atof(argv[++i]);
			c1 = atof(argv[++i]);
			c2 = atof(argv[++i]);
			c3 = atof(argv[++i]);
			break;
		case 'V':
			calc_vill = 1;
			break;
		case 'G':
			cut_view = 1;
			cut_view_type= atof(argv[++i]);
 			break;
		case 'g':
			cut_view = 2;
			cut_view_type= atof(argv[++i]);
 			break;

			/*case 'v':
			printf("evalglare  %s \n",version);
			exit(1); */

		case '1':
			output = 1;
			break;

		case 'v':
			if (argv[i][2] == '\0') {
                             printf("%s \n",RELEASENAME);				
			     exit(1);
			}
			if (argv[i][2] != 'f')
				goto userr;
			rval = viewfile(argv[++i], &userview, NULL);
			if (rval < 0) {
				fprintf(stderr,
						"%s: cannot open view file \"%s\"\n",
						progname, argv[i]);
				exit(1);
			} else if (rval == 0) {
				fprintf(stderr,
						"%s: bad view file \"%s\"\n", progname, argv[i]);
				exit(1);
			} else {
				gotuserview++;
			}
			break;


		default:
			goto userr;
		}
	}

/*fast calculation, if gendgp_profile is used: No Vertical illuminance calculation, only dgp is calculated*/

if (output == 1 && ext_vill == 1) {
                       calcfast=1;
		       }

/* read picture file */
	if (i == argc) {
		SET_FILE_BINARY(stdin);
		FILE *fp = fdopen(fileno(stdin), "rb");
		if (!(fp)) {
			fprintf(stderr, "fdopen on stdin failed\n");
			return EXIT_FAILURE;
		}
		if (!(pict_read_fp(p, fp)))
			return EXIT_FAILURE;;
	} else {
		if (!pict_read(p, argv[i]))
			return EXIT_FAILURE;
	}
	if (gotuserview) {
		if (p->valid_view)
			fprintf(stderr,
					"warning: overriding image view by commandline argument\n");
		if ((userview.horiz == 0) || (userview.vert == 0) || (userview.type == 0)) {
			fprintf(stderr, "error: if specified, a view must at least contain -vt -vv and -vh\n");
			return EXIT_FAILURE;
		}
		p->view = userview;
		if (!(pict_update_view(p))) {
			fprintf(stderr, "error: invalid view specified");
			return EXIT_FAILURE;
		}
		pict_update_evalglare_caches(p);
	} else if (!(p->valid_view)) {
		fprintf(stderr, "error: no valid view specified\n");
		return EXIT_FAILURE;
	}

/*		  fprintf(stderr,"Picture read!");*/

/* command line for output picture*/

	pict_set_comment(p, cline);

	if (p->view.type == VT_PAR) {
		fprintf(stderr, "wrong view type! must not be parallel ! \n");
		exit(1);
	}


/* Check size of search radius */
	rmx = (pict_get_xsize(p) / 2);
	rmy = (pict_get_ysize(p) / 2);
	rx = (pict_get_xsize(p) / 2 + 10);
	ry = (pict_get_ysize(p) / 2 + 10);
	r_center =
		acos(DOT(pict_get_cached_dir(p, rmx, rmy),
				 pict_get_cached_dir(p, rx, ry))) * 2 / 10;
	search_pix = max_angle / r_center;
	if (search_pix < 1.0) {
		fprintf(stderr,
				"warning: search radius less than 1 pixel! deactivating smoothing and peak extraction...\n");
		splithigh = 0;
		sgs = 0;

	} else {
		if (search_pix < 3.0) {
			fprintf(stderr,
					"warning: search radius less than 3 pixels! -> %f \n",
					search_pix);

		}
	}

/* Check task position  */

	if (task_lum == 1) {
		if (xt >= pict_get_xsize(p) || yt >= pict_get_ysize(p) || xt < 0
			|| yt < 0) {
			fprintf(stderr,
					"error: task position outside picture!! exit...");
			exit(1);
		}
	}


/*	printf("resolution: %dx%d\n",pict_get_xsize(p),pict_get_ysize(p)); */

	sang = 0.0;
	E_v = 0.0;
	E_v_dir = 0.0;
	avlum = 0.0;
	pict_new_gli(p);
	igs = 0;

/* cut out GUTH field of view and exit without glare evaluation */
if (cut_view==2) {
	if (cut_view_type==1) {
                   cut_view_1(p);
		   }else{
		   
		   	if (cut_view_type==2) {
                   	cut_view_2(p);
		   				}else{
		   					if (cut_view_type==3) {
		   					fprintf(stderr,"warning: pixel luminance is weighted by position index - do not use image for glare evaluations!!");
                   					cut_view_3(p);
            		   							}else{
		   					fprintf(stderr,"error: no valid option for view cutting!!");
                					exit(1); 
			   }
                   }}
		pict_write(p, file_out);
                exit(1);          }	






/* write positionindex into checkfile and exit, activated by option -p */
	if (posindex_picture == 1) {
		for (x = 0; x < pict_get_xsize(p); x++)
			for (y = 0; y < pict_get_ysize(p); y++) {
				if (pict_get_hangle
					(p, x, y, p->view.vdir, p->view.vup, &ang)) {
					if (pict_is_validpixel(p, x, y)) {
						lum =
							get_posindex(p, x, y,
										 posindex_2) / WHTEFFICACY;

						pict_get_color(p, x, y)[RED] = lum;
						pict_get_color(p, x, y)[GRN] = lum;
						pict_get_color(p, x, y)[BLU] = lum;
						lum_task = luminance(pict_get_color(p, x, y));
/*printf("x,y,posindex=%i %i %f %f\n",x,y,lum*WHTEFFICACY,lum_task);*/
					}
				}
			}
		pict_write(p, file_out);
		exit(1);
	}


/* calculation of vertical illuminance, average luminance, in case of filling activated-> fill picture */
/* fill, if necessary from 0 to yfillmin */

	if (fill == 1) {
		for (x = 0; x < pict_get_xsize(p); x++)
			for (y = yfillmin; y > 0; y = y - 1) {
				y1 = y + 1;
				lum = luminance(pict_get_color(p, x, y1));
				pict_get_color(p, x, y)[RED] = lum / 179.0;
				pict_get_color(p, x, y)[GRN] = lum / 179.0;
				pict_get_color(p, x, y)[BLU] = lum / 179.0;
			}
	}

	if (calcfast == 0) {
	for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
					lum = luminance(pict_get_color(p, x, y));
					if (fill == 1 && y >= yfillmax) {
						y1 = y - 1;
						lum = luminance(pict_get_color(p, x, y1));
						pict_get_color(p, x, y)[RED] = lum / 179.0;
						pict_get_color(p, x, y)[GRN] = lum / 179.0;
						pict_get_color(p, x, y)[BLU] = lum / 179.0;
					}

					if (lum > abs_max) {
						abs_max = lum;
					}
/* set luminance restriction, if -m is set */
					if (set_lum_max == 1 && lum >= lum_max) {
						pict_get_color(p, x, y)[RED] = new_lum_max / 179.0;
						pict_get_color(p, x, y)[GRN] = new_lum_max / 179.0;
						pict_get_color(p, x, y)[BLU] = new_lum_max / 179.0;
/*                                    printf("old lum, new lum %f %f \n",lum,luminance(pict_get_color(p,x,y))); */
						lum = luminance(pict_get_color(p, x, y));

					}
					if (set_lum_max2 == 1 && lum >= lum_max) {

						E_v_contr +=
							DOT(p->view.vdir,
								pict_get_cached_dir(p, x,
													y)) * pict_get_omega(p,
																		 x,
																		 y)
							* lum;
						omega_cos_contr +=
							DOT(p->view.vdir,
								pict_get_cached_dir(p, x,
													y)) * pict_get_omega(p,
																		 x,
																		 y)
							* 1;
					}


					sang += pict_get_omega(p, x, y);
					E_v +=
						DOT(p->view.vdir,
							pict_get_cached_dir(p, x,
												y)) * pict_get_omega(p, x,
																	 y) *
						lum;
					avlum +=
						pict_get_omega(p, x,
									   y) * luminance(pict_get_color(p, x,
																	 y));
				} else {
					pict_get_color(p, x, y)[RED] = 0.0;
					pict_get_color(p, x, y)[GRN] = 0.0;
					pict_get_color(p, x, y)[BLU] = 0.0;

				}
			}
		}

	if (set_lum_max2 == 1 && E_v_contr > 0 && (E_vl_ext - E_v) > 0) {

		lum_ideal = (E_vl_ext - E_v + E_v_contr) / omega_cos_contr;
		if (lum_ideal > 0 && lum_ideal < setvalue) {
			setvalue = lum_ideal;
		}
		printf
			("change luminance values!! lum_ideal,setvalue,E_vl_ext,E_v,E_v_contr %f  %f %f %f %f\n",
			 lum_ideal, setvalue, E_vl_ext, E_v, E_v_contr);


		for (x = 0; x < pict_get_xsize(p); x++)
			for (y = 0; y < pict_get_ysize(p); y++) {
				if (pict_get_hangle
					(p, x, y, p->view.vdir, p->view.vup, &ang)) {
					if (pict_is_validpixel(p, x, y)) {
						lum = luminance(pict_get_color(p, x, y));
						if (set_lum_max2 == 1 && lum >= lum_max) {

							pict_get_color(p, x, y)[RED] =
								setvalue / 179.0;
							pict_get_color(p, x, y)[GRN] =
								setvalue / 179.0;
							pict_get_color(p, x, y)[BLU] =
								setvalue / 179.0;

						}
					}
				}
			}

		pict_write(p, file_out2);

	}
	if (set_lum_max == 1) {
		pict_write(p, file_out2);

	}

	if (calc_vill == 1) {
		printf("%f\n", E_v);
		exit(1);
	}
	}else{
	/* in fast calculation mode: ev=ev_ext and sang=2*pi */
	sang=2*3.14159265359;
	lum_task =E_vl_ext/sang;
	E_v=E_vl_ext;
/*		printf("calc fast!! %f %f\n", sang,lum_task);*/
	
	
	}

/* cut out GUTH field of view for glare evaluation */
if (cut_view==1) {
	if (cut_view_type==1) {
                   cut_view_1(p);
		   }else{
		   
		   	if (cut_view_type==2) {
                   	cut_view_2(p);
		   				}else{
		   					if (cut_view_type==3) {
		   					fprintf(stderr,"warning: pixel luminance is weighted by position index - do not use image for glare evaluations!!");
                   					cut_view_3(p);
            		   							}else{
		   					fprintf(stderr,"error: no valid option for view cutting!!");
                					exit(1); 
			   }
                   }}
          }	

	if (calcfast == 0) {
	avlum = avlum / sang;
	lum_task = avlum;
	}
/*        if (ext_vill==1) {
		     E_v=E_vl_ext;
	             avlum=E_v/3.1415927;
                        }      */

	if (task_lum == 1) {
		lum_task = get_task_lum(p, xt, yt, omegat, task_color);
	}
	lum_source = lum_thres * lum_task;
/*         printf("Task Luminance=%f\n",lum_task);
                   pict_write(p,"task.pic");*/

	if (lum_thres > 100) {
		lum_source = lum_thres;
	}

	/*     printf("Ev, avlum, Omega_tot%f %f %f \n",E_v, avlum, sang); */

/* first glare source scan: find primary glare sources */
	for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
					act_lum = luminance(pict_get_color(p, x, y));
					if (act_lum > lum_source) {
						if (act_lum >lum_total_max) {
						                             lum_total_max=act_lum;
						                               }
						
						actual_igs =
							find_near_pgs(p, x, y, max_angle, 0, igs, 1);
						if (actual_igs == 0) {
							igs = igs + 1;
							pict_new_gli(p);
							pict_get_lum_min(p, igs) = HUGE_VAL;
							pict_get_Eglare(p,igs) = 0.0;
							pict_get_Dglare(p,igs) = 0.0;
							actual_igs = igs;
						}
						icol = setglcolor(p, x, y, actual_igs, uniform_gs, u_r, u_g , u_b);
						pict_get_gsn(p, x, y) = actual_igs;
						pict_get_pgs(p, x, y) = 1;
						add_pixel_to_gs(p, x, y, actual_igs);

					} else {
						pict_get_pgs(p, x, y) = 0;
					}
				}
			}
		}
/* 	pict_write(p,"firstscan.pic");   */

if (calcfast == 1 ) {
   skip_second_scan=1;
   }

/* second glare source scan: combine glare sources facing each other */
	change = 1;
	while (change == 1 && skip_second_scan==0) {
		change = 0;
		for (x = 0; x < pict_get_xsize(p); x++)
			for (y = 0; y < pict_get_ysize(p); y++) {
				before_igs = pict_get_gsn(p, x, y);
				if (pict_get_hangle
					(p, x, y, p->view.vdir, p->view.vup, &ang)) {
					if (pict_is_validpixel(p, x, y) && before_igs > 0) {
						actual_igs =
							find_near_pgs(p, x, y, max_angle, before_igs,
										  igs, 1);
						if (!(actual_igs == before_igs)) {
							change = 1;
						}
						if (before_igs > 0) {
							actual_igs = pict_get_gsn(p, x, y);
							icol = setglcolor(p, x, y, actual_igs, uniform_gs, u_r, u_g , u_b);
						}

					}
				}
			}
/*	pict_write(p,"secondscan.pic");*/

	}

/*smoothing: add secondary glare sources */
	i_max = igs;
	if (sgs == 1) {

/* simplified search radius: constant for entire picture, alway circle not an angle due to performance */
		x = (pict_get_xsize(p) / 2);
		y = (pict_get_ysize(p) / 2);
		rx = (pict_get_xsize(p) / 2 + 10);
		ry = (pict_get_ysize(p) / 2 + 10);
		r_center =
			acos(DOT(pict_get_cached_dir(p, x, y),
					 pict_get_cached_dir(p, rx, ry))) * 2 / 10;
		search_pix = max_angle / r_center / 1.75;

		for (x = 0; x < pict_get_xsize(p); x++) {
			for (y = 0; y < pict_get_ysize(p); y++) {
				if (pict_get_hangle
					(p, x, y, p->view.vdir, p->view.vup, &ang)) {
					if (pict_is_validpixel(p, x, y)
						&& pict_get_gsn(p, x, y) == 0) {
						for (i = 1; i <= i_max; i++) {

							if (pict_get_npix(p, i) > 0) {
								add_secondary_gs(p, x, y, search_pix, i, uniform_gs, u_r, u_g , u_b);
							}
						}

					}
				}
			}
		}

	}

/* extract extremes from glare sources to extra glare source */
	if (splithigh == 1 && lum_total_max>limit) {

		r_split = max_angle / 2.0;
		for (i = 0; i <= i_max; i++) {

			if (pict_get_npix(p, i) > 0) {
				l_max = pict_get_lum_max(p, i);
				i_splitstart = igs + 1;
				if (l_max >= limit) {
					for (x = 0; x < pict_get_xsize(p); x++)
						for (y = 0; y < pict_get_ysize(p); y++) {
							if (pict_get_hangle
								(p, x, y, p->view.vdir, p->view.vup, &ang))
							{
								if (pict_is_validpixel(p, x, y)
									&& luminance(pict_get_color(p, x, y))
									>= limit
									&& pict_get_gsn(p, x, y) == i) {
									if (i_splitstart == (igs + 1)) {
										igs = igs + 1;
										pict_new_gli(p);
										pict_get_Eglare(p,igs) = 0.0;
										pict_get_Dglare(p,igs) = 0.0;
										pict_get_lum_min(p, igs) =
											99999999999999.999;
										i_split = igs;
									} else {
										i_split =
											find_split(p, x, y, r_split,
													   i_splitstart, igs);
									}
									if (i_split == 0) {
										igs = igs + 1;
										pict_new_gli(p);
										pict_get_Eglare(p,igs) = 0.0;
										pict_get_Dglare(p,igs) = 0.0;
										pict_get_lum_min(p, igs) =
											99999999999999.999;
										i_split = igs;
									}
									split_pixel_from_gs(p, x, y, i_split, uniform_gs, u_r, u_g , u_b);

								}
							}
						}

				}
				change = 1;
				while (change == 1) {
					change = 0;
					for (x = 0; x < pict_get_xsize(p); x++)
						for (y = 0; y < pict_get_ysize(p); y++) {
							before_igs = pict_get_gsn(p, x, y);
							if (before_igs >= i_splitstart) {
								if (pict_get_hangle
									(p, x, y, p->view.vdir, p->view.vup,
									 &ang)) {
									if (pict_is_validpixel(p, x, y)
										&& before_igs > 0) {
										actual_igs =
											find_near_pgs(p, x, y,
														  max_angle,
														  before_igs, igs,
														  i_splitstart);
										if (!(actual_igs == before_igs)) {
											change = 1;
										}
										if (before_igs > 0) {
											actual_igs =
												pict_get_gsn(p, x, y);
											icol =
												setglcolor(p, x, y,
														   actual_igs, uniform_gs, u_r, u_g , u_b);
										}

									}
								}
							}

						}
				}


			}
		}
	}

/* calculation of direct vertical illuminance for CGI and for disability glare*/


	if (calcfast == 0) {
	for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
					if (pict_get_gsn(p, x, y) > 0) {
						pict_get_Eglare(p,pict_get_gsn(p, x, y)) +=
							DOT(p->view.vdir, pict_get_cached_dir(p, x, y))
							* pict_get_omega(p, x, y)
							* luminance(pict_get_color(p, x, y));
						E_v_dir +=
							DOT(p->view.vdir, pict_get_cached_dir(p, x, y))
							* pict_get_omega(p, x, y)
							* luminance(pict_get_color(p, x, y));
					}
				}
			}
		}
	lum_backg_cos = (E_v - E_v_dir) / 3.1415927;
	lum_backg = lum_backg_cos;
         }
/* calc of band luminance if applied */
	if (band_calc == 1) {
        band_avlum=get_band_lum( p, band_angle, band_color);
	}

/*printf("total number of glare sources: %i\n",igs); */
	lum_sources = 0;
	omega_sources = 0;
	for (x = 0; x <= igs; x++) {
		lum_sources += pict_get_av_lum(p, x) * pict_get_av_omega(p, x);
		omega_sources += pict_get_av_omega(p, x);
	}
	if (non_cos_lb == 1) {
		lum_backg =
			(sang * avlum - lum_sources) / (sang - omega_sources);
	}
/* print detailed output */
	if (detail_out == 1) {
		i = 0;
		for (x = 0; x <= igs; x++) {
/* resorting glare source numbers */
			if (pict_get_npix(p, x) > 0) {
				i = i + 1;
			}
		}



		printf
			("%i No pixels x-pos y-pos L_s Omega_s Posindx L_b L_t E_vert Edir Max_Lum Sigma xdir ydir zdir Eglare_cie Lveil_cie \n",
			 i);
		if (i == 0) {
			printf("%i %f %f %f %f %.10f %f %f %f %f %f %f\n", i, 0.0, 0.0,
				   0.0, 0.0, 0.0, 0.0, lum_backg, lum_task, E_v, E_v_dir,
				   abs_max);

		} else {
			i = 0;
			for (x = 0; x <= igs; x++) {
				if (pict_get_npix(p, x) > 0) {
					i = i + 1;
					pict_get_sigma(p, pict_get_av_posx(p, x),
								   pict_get_av_posy(p, x), p->view.vdir,
								   p->view.vup, &sigma);
								   
					x2=pict_get_av_posx(p, x);
					y2=pict_get_av_posy(p, x);
					teta = 180.0 / 3.1415927 * acos(DOT(p->view.vdir, pict_get_cached_dir(p, x2, y2))); 
					Lveil_cie = 10*pict_get_Eglare(p,x)/(teta*teta+0.0000000000000001);

	if (apply_disability == 1 && Lveil_cie <=disability_thresh) { 
	                                                             Lveil_cie =0 ;
	                                                           }
					Lveil_cie_sum =  Lveil_cie_sum + Lveil_cie; 
					printf("%i %f %f %f %f %.10f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
						   i, pict_get_npix(p, x), pict_get_av_posx(p, x),
						   pict_get_ysize(p) - pict_get_av_posy(p, x),
						   pict_get_av_lum(p, x), pict_get_av_omega(p, x),
						   get_posindex(p, pict_get_av_posx(p, x),
										pict_get_av_posy(p, x),
										posindex_2), lum_backg, lum_task,
						   E_v, E_v_dir, abs_max, sigma * 180 / 3.1415927
						   ,pict_get_cached_dir(p, x2,y2)[0],pict_get_cached_dir(p, x2,y2)[1],pict_get_cached_dir(p, x2,y2)[2],pict_get_Eglare(p,x),Lveil_cie,teta
						   
						   );
				}
			}
		}
	}



/* calculation of indicees */

/* check vertical illuminance range */
	E_v2 = E_v;

	if (E_v2 < 100) {
		fprintf(stderr,
				"Notice: Vertical illuminance is below 100 lux !!\n");
	}
	dgp =
		get_dgp(p, E_v2, igs, a1, a2, a3, a4, a5, c1, c2, c3, posindex_2);
/* low light correction */
       if (E_v < 1000) {
       low_light_corr=1.0*exp(0.024*E_v-4)/(1+exp(0.024*E_v-4));} else {low_light_corr=1.0 ;}
       dgp =low_light_corr*dgp;
       
/* age correction */
       
	if (age_corr == 1) {
	         age_corr_factor=1.0/(1.1-0.5*age/100.0);
		 }
       dgp =age_corr_factor*dgp;
       
       
	if (ext_vill == 1) {
		if (E_vl_ext < 100) {
			fprintf(stderr,
					"Notice: Vertical illuminance is below 100 lux !!\n");
		}
	}

if (calcfast == 0) {
	dgi = get_dgi(p, lum_backg, igs, posindex_2);
	ugr = get_ugr(p, lum_backg, igs, posindex_2);
	cgi = get_cgi(p, E_v, E_v_dir, igs, posindex_2);
	vcp = get_vcp(p, avlum, igs, posindex_2);
	Lveil = get_disability(p, avlum, igs);
}
/* check dgp range */
	if (dgp <= 0.2) {
		fprintf(stderr,
				"Notice: Low brightness scene. dgp below 0.2! dgp might underestimate glare sources\n");
	}
	if (E_v < 380) {
		fprintf(stderr,
				"Notice: Low brightness scene. Vertical illuminance less than 380 lux! dgp might underestimate glare sources\n");
	}



	if (output == 0) {
		if (detail_out == 1) {
			if (ext_vill == 1) {
	        dgp_ext =
		               get_dgp(p, E_vl_ext, igs, a1, a2, a3, a4, a5, c1, c2, c3,
				posindex_2);
				dgp = dgp_ext;
				if (E_vl_ext < 1000) {
				low_light_corr=1.0*exp(0.024*E_vl_ext-4)/(1+exp(0.024*E_vl_ext-4));} else {low_light_corr=1.0 ;}
				dgp =low_light_corr*dgp;
       				dgp =age_corr_factor*dgp;
			}

			printf
				("dgp,av_lum,E_v,lum_backg,E_v_dir,dgi,ugr,vcp,cgi,lum_sources,omega_sources,Lveil,Lveil_cie,band_avlum: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",
				 dgp, avlum, E_v, lum_backg, E_v_dir, dgi, ugr, vcp, cgi,
				 lum_sources / omega_sources, omega_sources, Lveil,Lveil_cie_sum,band_avlum);
		} else {
			if (detail_out2 == 1) {

				printf
					("dgp,dgi,ugr,vcp,cgi,dgp_ext,Ev_calc,abs_max,Lveil: %f %f %f %f %f %f %f %f %f \n",
					 dgp, dgi, ugr, vcp, cgi, dgp_ext, E_v, abs_max,
					 Lveil);

			} else {
				if (ext_vill == 1) {
	        dgp_ext = get_dgp(p, E_vl_ext, igs, a1, a2, a3, a4, a5, c1, c2, c3,posindex_2);
				
				if (E_vl_ext < 1000) {
				low_light_corr=1.0*exp(0.024*E_vl_ext-4)/(1+exp(0.024*E_vl_ext-4));} else {low_light_corr=1.0 ;} 
					dgp =low_light_corr*dgp_ext;
       					dgp =age_corr_factor*dgp;
				}
				printf("dgp,dgi,ugr,vcp,cgi,Lveil: %f %f %f %f %f %f  \n",
					   dgp, dgi, ugr, vcp, cgi, Lveil);

			}
		}
	} else {
	        dgp_ext =
		               get_dgp(p, E_vl_ext, igs, a1, a2, a3, a4, a5, c1, c2, c3,
				posindex_2);
				dgp = dgp_ext;
				if (E_vl_ext < 1000) {
				low_light_corr=1.0*exp(0.024*E_vl_ext-4)/(1+exp(0.024*E_vl_ext-4)); } else {low_light_corr=1.0 ;} 
				dgp =low_light_corr*dgp;
       				dgp =age_corr_factor*dgp;
		printf("%f\n", dgp);
	}


/*printf("lowlight=%f\n", low_light_corr); */


/*	printf("hallo \n");


			apply_disability = 1;
			disability_thresh = atof(argv[++i]);
coloring of disability glare sources red, remove all other colors



this function was removed because of bugs.... 
has to be re-written from scratch....


*/








/*output  */
	if (checkfile == 1) {
		pict_write(p, file_out);
	}



	return EXIT_SUCCESS;
	exit(0);

  userr:
	fprintf(stderr,
			"Usage: %s [-s][-d][-c picture][-t xpos ypos angle] [-T xpos ypos angle] [-b fact] [-r angle] [-y] [-Y lum] [-i Ev] [-I Ev ymax ymin] [-v] picfile\n",
			progname);
	exit(1);
}



#endif
