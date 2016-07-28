#ifndef lint
static const char RCSid[] = "$Id: evalglare.c,v 2.4 2016/07/28 16:46:10 greg Exp $";
#endif
/* EVALGLARE V1.29
 * Evalglare Software License, Version 2.0
 *
 * Copyright (c) 1995 - 2016 Fraunhofer ISE, EPFL.
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
                developed at Fraunhofer ISE and EPFL by J. Wienold" and 
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Evalglare," "EPFL" and "Fraunhofer ISE" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact jan.wienold@epfl.ch
 *
 * 5. Products derived from this software may not be called "evalglare",
 *       nor may "evalglare" appear in their name, without prior written
 *       permission of EPFL.
 *
 * Redistribution and use in source and binary forms, WITH 
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *
 * conditions 1.-5. apply
 * 
 * 6. In order to ensure scientific correctness, any modification in source code imply fulfilling all following comditions:
 *     a) A written permission from EPFL. For written permission, please contact jan.wienold@epfl.ch.
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

/* evalglare.c, v1.18 2015/11/10 add masking, UDP, PSGV, DGI_mod, UGR_exp and zoning. For zoning and band luminance, a statistical package from C. Reetz was added in order to derive median, std_dev, percentiles.
memory leakage was checked and is o.k.
   */

/* evalglare.c, v1.19 2015/12/04 sorting subroutines in statistical package changed by. C. Reetz, adding Position index weighted Luminance evaluation (mean and median of entire image, only in detailed output available)
   */
/* evalglare.c, v1.20 2015/12/21 removing -ansi gcc-option in makefile for the standalone version, adding DGR as output since VCP equation seems not to match with original publication
   */
/* evalglare.c, v1.21 2016/03/15 add a new pixel-overflow-correction option -N 
   */
/* evalglare.c, v1.22 2016/03/22 fixed problem using -N option and angular distance calculation, when actual pixel equals disk-center (arccos returned nan)
   */
/* evalglare.c, v1.23 2016/04/18 fixed problem making black corners for -vth fish-eye images
   */
/* evalglare.c, v1.24 2016/04/26 significant speed improvement for 2nd glare source scan: for other glare sources will not be searched any more, when the 8 surrounding pixels have the same glare source number! Smaller speed improvement in the first glare source scan: remembering if pixel before was also a glare source, then no search for surrounding glare sources is applied.
changed masking threshold to 0.05 cd/m2
   */
/* evalglare.c, v1.25 2016/04/27 removed the change on the first glare source scan: in few cases glare sources are not merged together in the same way as before. Although the difference is marginal, this algorithm was removed in order to be consistent to the existing results. 
   */
/* evalglare.c, v1.26 2016/04/28 removed the bug Lb getting nan in case all glare source pixels are above the peak extraction limit. 
 The calculation of the average source lumiance and solid angle was adapted to be more robust for extreme cases.
 In case no glare source is found, the values of the glare metrics, relying only on glare sources, is set to zero (100 for VCP) to avoid nan and -inf in the output. 
 Changed glare section output when 0 glare source are found to have the same amount of columns than when glare sources are found
   */
/* evalglare.c, v1.27 2016/06/13 include a warning, if -vtv is in the header. Check, if the corners are black AND -vtv is used ->stopping (stopping can be avoided by using the forcing option -f ). Check, if an invalid exposure is in the header -> stopping in any case. 
   */
/* evalglare.c, v1.28 2016/07/08 small code changes in the section of calculating the E_glare (disability) for each glare source. 
   */
/* evalglare.c, v1.29 2016/07/12 change threshold for the black corner to 2cd/m2. 
   */
/* evalglare.c, v1.30 2016/07/28 change zonal output: If just one zone is specified, only one zone shows up in the output (bug removal). 
   */
#define EVALGLARE
#define PROGNAME "evalglare"
#define VERSION "1.30 release 29.07.2016 by EPFL, J.Wienold"
#define RELEASENAME PROGNAME " " VERSION


#include "pictool.h"
#include "rtio.h"
#include <math.h>
#include <string.h>
#include "platform.h"
#include "muc_randvar.h"

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





/* subroutine for coloring the glare sources 
     red is used only for explicit call of the subroutine with acol=0  , e.g. for disability glare
     -1: set to grey*/
int setglcolor(pict * p, int x, int y, int acol, int uniform_gs, double u_r, double u_g ,double u_b)
{
	int icol;
	double pcol[3][1000], act_lum, l;

        l=u_r+u_g+u_b ;
	
	pcol[RED][0] = 1.0 / CIE_rf;
	pcol[GRN][0] = 0.0 / CIE_gf;
	pcol[BLU][0] = 0.0 / CIE_bf;

	pcol[RED][1] = 0.0 / CIE_rf;
	pcol[GRN][1] = 1.0 / CIE_gf;
	pcol[BLU][1] = 0.0 / CIE_bf;

	pcol[RED][2] = 0.15 / CIE_rf;
	pcol[GRN][2] = 0.15 / CIE_gf;
	pcol[BLU][2] = 0.7 / CIE_bf;

	pcol[RED][3] = .5 / CIE_rf;
	pcol[GRN][3] = .5 / CIE_gf;
	pcol[BLU][3] = 0.0 / CIE_bf;

	pcol[RED][4] = .5 / CIE_rf;
	pcol[GRN][4] = .0 / CIE_gf;
	pcol[BLU][4] = .5 / CIE_bf ;

	pcol[RED][5] = .0 / CIE_rf;
	pcol[GRN][5] = .5 / CIE_gf;
	pcol[BLU][5] = .5 / CIE_bf;

	pcol[RED][6] = 0.333 / CIE_rf;
	pcol[GRN][6] = 0.333 / CIE_gf;
	pcol[BLU][6] = 0.333 / CIE_bf;


	pcol[RED][999] = 1.0 / WHTEFFICACY;
	pcol[GRN][999] = 1.0 / WHTEFFICACY;
	pcol[BLU][999] = 1.0 / WHTEFFICACY;
 
 
 	pcol[RED][998] = u_r /(l* CIE_rf) ;
	pcol[GRN][998] = u_g /(l* CIE_gf);
	pcol[BLU][998] = u_b /(l* CIE_bf);
/*        printf("CIE_rf,CIE_gf,CIE_bf,l=%f %f %f %f %f\n",CIE_rf,CIE_gf,CIE_bf,l,WHTEFFICACY); */
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
/*        printf("x,y,lum_before,lum_after,diff %i %i %f %f %f\n",x,y, act_lum,luminance(pict_get_color(p, x, y)), act_lum-luminance(pict_get_color(p, x, y))); */
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

/*		icol = setglcolor(p, x, y, act_gs, uniform_gs, u_r, u_g , u_b); */

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

/*	icol = setglcolor(p, x, y, new_gsn, uniform_gs, u_r, u_g , u_b); */
        

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

        if (postype == 1) {
/* KIM  model */        
        posindex = exp ((sigma-(-0.000009*tau*tau*tau+0.0014*tau*tau+0.0866*tau+21.633))/(-0.000009*tau*tau*tau+0.0013*tau*tau+0.0853*tau+8.772));
        }else{
/* Guth model, equation from IES lighting handbook */
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
}

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

/* subroutine for the calculation of the daylight glare index DGI*/


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
/* subroutine for the calculation of the modified daylight glare index DGI_mod (Fisekis2003) 
   several glare sources are taken into account and summed up, original equation has no summary
*/


float get_dgi_mod(pict * p, float lum_a, int igs, int posindex_2)
{
	float dgi_mod, sum_glare, omega_s;
	int i;


	sum_glare = 0;
	omega_s = 0;

	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			omega_s =
				pict_get_av_omega(p, i) / get_posindex(p,pict_get_av_posx(p,i), pict_get_av_posy(p,i),posindex_2) /
				get_posindex(p, pict_get_av_posx(p, i), pict_get_av_posy(p, i), posindex_2);
			sum_glare += 0.478 * pow(pict_get_av_lum(p, i), 1.6) * pow(omega_s,0.8) / (pow(lum_a,0.85) + 0.07 * pow(pict_get_av_omega(p, i),0.5) * pict_get_av_lum(p, i));
                    /*    printf("i,sum_glare %i %f\n",i,sum_glare); */
		}
	}
	dgi_mod = 10 * log10(sum_glare);

	return dgi_mod;

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
float get_dgr(pict * p, double lum_a, int igs, int posindex_2)
{
	float dgr;
	double sum_glare;
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



	return dgr;

}

float get_vcp(float dgr )
{
	float vcp;

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
/* subroutine for the calculation of the experimental unified glare rating (Fisekis 2003)*/
float get_ugr_exp(pict * p, double lum_backg, double lum_a, int igs, int posindex_2)
{
	float ugr_exp;
	double sum_glare;
	int i, i_glare;


	sum_glare = 0;
	i_glare = 0;
	for (i = 0; i <= igs; i++) {

		if (pict_get_npix(p, i) > 0) {
			i_glare = i_glare + 1;
			sum_glare +=
				pow(1 /get_posindex(p, pict_get_av_posx(p, i), pict_get_av_posy(p, i), posindex_2),2) *pict_get_av_lum(p, i)* pict_get_av_omega(p, i);

		}
	}
	ugr_exp =8 * log10(lum_a) + 8 * log10(sum_glare/lum_backg);

	return ugr_exp;

}
/* subroutine for the calculation of the unified glare probability (Hirning 2014)*/
float get_ugp(pict * p, double lum_backg, int igs, int posindex_2)
{
	float ugp;
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
	ugp = 0.26 * log10(0.25 / lum_backg * sum_glare);

	return ugp;

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

float get_cgi(pict * p, double E_v, double E_v_dir, int igs, int posindex_2)
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

/* subroutine for the calculation of the PGSV; is only applied, when masking is done, since it refers only to the window. Important: masking area must be the window! */
float get_pgsv(double E_v, double E_mask,double omega_mask,double lum_mask_av)
{
	float pgsv;
	double Lb;

        Lb = (E_v-E_mask)/1.414213562373;

	pgsv =3.2*log10(lum_mask_av)-0.64*log10(omega_mask)+(0.79*log10(omega_mask)-0.61)*log10(Lb)-8.2 ;


	return pgsv;
}

/* subroutine for the calculation of the PGSV_saturation; is only applied, when masking is done, since it refers only to the window. Important: masking area must be the window! */
float get_pgsv_sat(double E_v)
{
	float pgsv_sat;

	pgsv_sat =3.3-(0.57+3.3)/pow((1+E_v/1.414213562373/1250),1.7);


	return pgsv_sat;
}




#ifdef	EVALGLARE


/* main program 
------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char **argv)
{
	#define CLINEMAX 4095 /* memory allocated for command line string */
	pict *p = pict_create();
        pict *pm = pict_create();
	int     skip_second_scan,calcfast,age_corr,cut_view,cut_view_type,calc_vill, output, detail_out2, x1,y1, fill, yfillmax, yfillmin,
		ext_vill, set_lum_max, set_lum_max2, img_corr,x_disk,y_disk,task_color, i_splitstart,zones,act_gsn,splitgs,
		i_split, posindex_2, task_lum, checkfile, rval, i, i_max, x, y,x2,y2,x_zone,y_zone, i_z1, i_z2,
		igs, actual_igs, lastpixelwas_gs, icol, xt, yt, change,checkpixels, before_igs, sgs, splithigh,uniform_gs,x_max, y_max,y_mid,
		detail_out, posindex_picture, non_cos_lb, rx, ry, rmx,rmy,apply_disability,band_calc,band_color,masking,i_mask,no_glaresources,force;
	double  lum_total_max,age_corr_factor,age,dgp_ext,dgp,low_light_corr,omega_cos_contr, setvalue, lum_ideal, E_v_contr, sigma,om,delta_E,
		E_vl_ext, lum_max, new_lum_max, r_center, ugp, ugr_exp, dgi_mod,lum_a, pgsv,E_v_mask,pgsv_sat,angle_disk,dist,n_corner_px,zero_corner_px,
		search_pix, a1, a2, a3, a4, a5, c3, c1, c2, r_split, max_angle,r_actual,lum_actual,
		omegat, sang, E_v, E_v2, E_v_dir, avlum, act_lum, ang, angle_z1, angle_z2,per_95_band,per_75_band,pos,
		l_max, lum_backg, lum_backg_cos, omega_sources, lum_sources,per_75_mask,per_95_mask,per_75_z1,per_95_z1,per_75_z2,per_95_z2,
		lum, lum_source,teta,Lveil_cie,Lveil_cie_sum,disability_thresh,u_r,u_g,u_b,band_angle,band_avlum,
		lum_mask[1],lum_mask_av,lum_mask_std[1],lum_mask_median[1],omega_mask,bbox[2],
		lum_band[1],lum_band_av,lum_band_std[1],lum_band_median[1],omega_band,bbox_band[2],
		lum_z1[1],lum_z1_av,lum_z1_std[1],lum_z1_median[1],omega_z1,bbox_z1[2],
		lum_z2[1],lum_z2_av,lum_z2_std[1],lum_z2_median[1],omega_z2,bbox_z2[2],
		lum_pos[1],lum_nopos_median[1],lum_pos_median[1],lum_pos2_median[1],lum_pos_mean,lum_pos2_mean;
	float lum_task, lum_thres, dgi,  vcp, cgi, ugr, limit, dgr, 
		abs_max, Lveil;
	char maskfile[500],file_out[500], file_out2[500], version[500];
	char *cline;
	VIEW userview = STDVIEW;
	int gotuserview = 0;
	struct muc_rvar* s_mask;
	s_mask = muc_rvar_create();
        muc_rvar_set_dim(s_mask, 1);
	muc_rvar_clear(s_mask);
	struct muc_rvar* s_band;
	s_band = muc_rvar_create();
        muc_rvar_set_dim(s_band, 1);
	muc_rvar_clear(s_band);
	struct muc_rvar* s_z1;
	s_z1 = muc_rvar_create();
        muc_rvar_set_dim(s_z1, 1);
	muc_rvar_clear(s_z1);

	struct muc_rvar* s_z2;
	s_z2 = muc_rvar_create();
        muc_rvar_set_dim(s_z2, 1);
	muc_rvar_clear(s_z2);

	struct muc_rvar* s_noposweight;
	s_noposweight = muc_rvar_create();
        muc_rvar_set_dim(s_noposweight, 1);
	muc_rvar_clear(s_noposweight);
 
	struct muc_rvar* s_posweight;
	s_posweight = muc_rvar_create();
        muc_rvar_set_dim(s_posweight, 1);
	muc_rvar_clear(s_posweight);

	struct muc_rvar* s_posweight2;
	s_posweight2 = muc_rvar_create();
        muc_rvar_set_dim(s_posweight2, 1);
	muc_rvar_clear(s_posweight2);

	/*set required user view parameters to invalid values*/
        delta_E=0.0;
        no_glaresources=0;
        n_corner_px=0;
        zero_corner_px=0;
        force=0;
        dist=0.0;
        u_r=0.33333;
        u_g=0.33333;
        u_b=0.33333;
        band_angle=0;
        angle_z1=0;
        angle_z2=0;
        x_zone=0;
        y_zone=0;
        per_75_z2=0;
        per_95_z2=0;
        lum_pos_mean=0;
        lum_pos2_mean=0;
	lum_band_av = 0.0;
	omega_band = 0.0;
	pgsv = 0.0 ;
	pgsv_sat = 0.0 ;
	E_v_mask = 0.0;
	lum_z1_av = 0.0;
	omega_z1 = 0.0;
	lum_z2_av = 0.0;
	omega_z2 = 0.0 ;
	i_z1 = 0;
	i_z2 = 0;        
	zones = 0;
	lum_z2_av = 0.0;
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
	x_disk=0;
	y_disk=0;
	angle_disk=0;
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
	img_corr=0;
	abs_max = 0;
	progname = argv[0];
	E_v_contr = 0.0;
	strcpy(version, "1.19 release 09.12.2015 by J.Wienold");
        low_light_corr=1.0;
 	output = 0;
	calc_vill = 0;
	band_avlum = -99;
	band_calc = 0;
        masking = 0;
	lum_mask[0]=0.0;
	lum_mask_av=0.0;
	omega_mask=0.0;
	i_mask=0;
        actual_igs=0;
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
		case 'A':
			masking = 1;
			detail_out = 1;
			strcpy(maskfile, argv[++i]);
			pict_read(pm,maskfile );

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
		case 'f':
			force = 1;
			break;
		case 'k':
			apply_disability = 1;
			disability_thresh = atof(argv[++i]);
			break;
			
		case 'p':
			posindex_picture = 1;
			break;
		case 'P':
			posindex_2 = atoi(argv[++i]);
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
			img_corr = 1;
			set_lum_max = 1;
			lum_max = atof(argv[++i]);
			new_lum_max = atof(argv[++i]);
			strcpy(file_out2, argv[++i]);
/*			printf("max lum set to %f\n",new_lum_max);*/
			break;

		case 'M':
			img_corr = 1;
			set_lum_max2 = 1;
			lum_max = atof(argv[++i]);
			E_vl_ext = atof(argv[++i]);
			strcpy(file_out2, argv[++i]);
/*			printf("max lum set to %f\n",new_lum_max);*/
			break;
		case 'N':
			img_corr = 1;
			set_lum_max2 = 2;
			x_disk = atoi(argv[++i]);
			y_disk = atoi(argv[++i]);
			angle_disk = atof(argv[++i]);
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
		case 'l':
			zones = 1;
			x_zone = atoi(argv[++i]);
			y_zone = atoi(argv[++i]);
			angle_z1 = atof(argv[++i]);
			angle_z2 = -1;
			break;
		case 'L':
			zones = 2;
			x_zone = atoi(argv[++i]);
			y_zone = atoi(argv[++i]);
			angle_z1 = atof(argv[++i]);
			angle_z2 = atof(argv[++i]);
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
/*masking and zoning cannot be applied at the same time*/

if (masking ==1 && zones >0) {
	       fprintf(stderr,	" masking and zoning cannot be activated at the same time!\n");
               exit(1);
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

/* several checks */
	if (p->view.type == VT_PAR) {
		fprintf(stderr, "error: wrong view type! must not be parallel ! \n");
		exit(1);
	}

	if (p->view.type == VT_PER) {
		fprintf(stderr, "warning: image has perspective view type specified in header ! \n");
	}

	if (masking == 1) {

		if (!pict_get_xsize(p)==pict_get_xsize(pm) || !pict_get_ysize(p)==pict_get_ysize(pm)) {
		fprintf(stderr, "error: masking image has other resolution than main image ! \n");
		fprintf(stderr, "size must be identical \n");
	        printf("resolution main image : %dx%d\n",pict_get_xsize(p),pict_get_ysize(p)); 
	        printf("resolution masking image : %dx%d\n",pict_get_xsize(pm),pict_get_ysize(pm)); 
		exit(1);
		
		}

	}
/*	printf("resolution: %dx%d\n",pict_get_xsize(p),pict_get_ysize(p)); */

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
	pict_get_z1_gsn(p,igs) = 0;
	pict_get_z2_gsn(p,igs) = 0;


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
                   lum = luminance(pict_get_color(p, x, y));		  
                   dist=sqrt((x-rmx)*(x-rmx)+(y-rmy)*(y-rmy));
		  if (dist>((rmx+rmy)/2)) {
		     n_corner_px=n_corner_px+1;
		     if (lum<7.0) {
		         zero_corner_px=zero_corner_px+1;
		         }
		     }
		
		
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
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
					if (img_corr == 1 ) {
						if (set_lum_max == 1 && lum >= lum_max) {
							pict_get_color(p, x, y)[RED] = new_lum_max / 179.0;
							pict_get_color(p, x, y)[GRN] = new_lum_max / 179.0;
							pict_get_color(p, x, y)[BLU] = new_lum_max / 179.0;
/*                                    			printf("old lum, new lum %f %f \n",lum,luminance(pict_get_color(p,x,y))); */
							lum = luminance(pict_get_color(p, x, y));
							}
						if (set_lum_max2 == 1 && lum >= lum_max) {

							E_v_contr +=DOT(p->view.vdir,pict_get_cached_dir(p, x,y)) * pict_get_omega(p,x,y)* lum;
							omega_cos_contr +=DOT(p->view.vdir,pict_get_cached_dir(p, x,y)) * pict_get_omega(p,x,y)* 1;
							}
						if (set_lum_max2 == 2 ) {
							r_actual =acos(DOT(pict_get_cached_dir(p, x_disk, y_disk), pict_get_cached_dir(p, x, y))) * 2;
							if (x_disk == x && y_disk==y ) r_actual=0.0;

			                                act_lum = luminance(pict_get_color(p, x, y));

			                               if (r_actual <= angle_disk) {
							      E_v_contr +=DOT(p->view.vdir,pict_get_cached_dir(p, x,y)) * pict_get_omega(p,x,y)* lum;
				                              omega_cos_contr += DOT(p->view.vdir,pict_get_cached_dir(p, x,y)) * pict_get_omega(p,x,y)* 1; 
								}

						
						
							}
					}
					om = pict_get_omega(p, x, y);
					sang += om;
					E_v += DOT(p->view.vdir, pict_get_cached_dir(p, x,y)) * om *lum;
					avlum += om * lum;
					pos=get_posindex(p, x, y,posindex_2);

					lum_pos[0]=lum;
	                                muc_rvar_add_sample(s_noposweight,lum_pos);
					lum_pos[0]=lum/pos;
                                        lum_pos_mean+=lum_pos[0]*om;
	                                muc_rvar_add_sample(s_posweight, lum_pos);
					lum_pos[0]=lum_pos[0]/pos;
                                        lum_pos2_mean+=lum_pos[0]*om;
	                                muc_rvar_add_sample(s_posweight2,lum_pos);



				} else {
					pict_get_color(p, x, y)[RED] = 0.0;
					pict_get_color(p, x, y)[GRN] = 0.0;
					pict_get_color(p, x, y)[BLU] = 0.0;

				}
			}else {
					pict_get_color(p, x, y)[RED] = 0.0;
					pict_get_color(p, x, y)[GRN] = 0.0;
					pict_get_color(p, x, y)[BLU] = 0.0;

				}
		}

/* check if image has black corners AND a perspective view */

       if (zero_corner_px/n_corner_px > 0.70 && (p->view.type == VT_PER) ) {
       printf (" corner pixels are to  %f %% black! \n",zero_corner_px/n_corner_px*100);
 		printf("error! The image has black corners AND header contains a perspective view type definition !!!\n") ;
		
		if (force==0) {
				printf("stopping...!!!!\n") ;

                      exit(1);
                 }
       } 




 	lum_pos_mean= lum_pos_mean/sang;
 	lum_pos2_mean= lum_pos2_mean/sang;

	if (set_lum_max2 >= 1 && E_v_contr > 0 && (E_vl_ext - E_v) > 0) {

		lum_ideal = (E_vl_ext - E_v + E_v_contr) / omega_cos_contr;
		if (set_lum_max2 == 2 && lum_ideal >= 2e9) {
		printf("warning! luminance of replacement pixels would be larger than 2e9 cd/m2. Value set to 2e9cd/m2!\n") ;
		}

		if (lum_ideal > 0 && lum_ideal < setvalue) {
			setvalue = lum_ideal;
		}
		printf("change luminance values!! lum_ideal,setvalue,E_vl_ext,E_v,E_v_contr %f  %f %f %f %f\n",
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

						}else{ if(set_lum_max2 ==2 ) {
							r_actual =acos(DOT(pict_get_cached_dir(p, x_disk, y_disk), pict_get_cached_dir(p, x, y))) * 2;
							if (x_disk == x && y_disk==y ) r_actual=0.0;

			                               if (r_actual <= angle_disk) {
			                               
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
				}
			}
			

		pict_write(p, file_out2);
        exit(1);
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
	for (x = 0; x < pict_get_xsize(p); x++) { 
/*                lastpixelwas_gs=0; */
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
					act_lum = luminance(pict_get_color(p, x, y));
					if (act_lum > lum_source) {
						if (act_lum >lum_total_max) {
						                             lum_total_max=act_lum;
						                               }
/* speed improvement first scan: when last pixel was glare source, then use this glare source number instead of searching. 
   has been deactivated with v1.25
                       
						if (lastpixelwas_gs==0 || search_pix <= 1.0 ) { */
						actual_igs = find_near_pgs(p, x, y, max_angle, 0, igs, 1);
/*  }*/
						if (actual_igs == 0) {
							igs = igs + 1;
							pict_new_gli(p);
							pict_get_lum_min(p, igs) = HUGE_VAL;
							pict_get_Eglare(p,igs) = 0.0;
							pict_get_Dglare(p,igs) = 0.0;
							pict_get_z1_gsn(p,igs) = 0;
							pict_get_z2_gsn(p,igs) = 0;
							actual_igs = igs;
							
						}
/* no coloring any more here 		icol = setglcolor(p, x, y, actual_igs, uniform_gs, u_r, u_g , u_b); */
						pict_get_gsn(p, x, y) = actual_igs;
						pict_get_pgs(p, x, y) = 1;
						add_pixel_to_gs(p, x, y, actual_igs);
/*                                                lastpixelwas_gs=actual_igs; */

					} else {
						pict_get_pgs(p, x, y) = 0;
 /*                                               lastpixelwas_gs=0;*/
					}
				}
			}
		} 
             }
/* 	pict_write(p,"firstscan.pic");   */


if (calcfast == 1 || search_pix <= 1.0) {
   skip_second_scan=1;
   }

/* second glare source scan: combine glare sources facing each other */
	change = 1;
        i = 0;
	while (change == 1 && skip_second_scan==0) {
		change = 0;
                i = i+1;
		for (x = 0; x < pict_get_xsize(p); x++) {
			for (y = 0; y < pict_get_ysize(p); y++) {
				if (pict_get_hangle
					(p, x, y, p->view.vdir, p->view.vup, &ang)) {
                                        checkpixels=1;
  				        before_igs = pict_get_gsn(p, x, y);

/* speed improvement: search for other nearby glare sources only, when at least one adjacend pixel has another glare source number  */
   				        if (x >1 && x<pict_get_xsize(p)-2 && y >1 && y<pict_get_ysize(p)-2 )   {
                                                x1=x-1;
                                                x2=x+1;
                                                y1=y-1;
                                                y2=y+1;
                                            if (before_igs == pict_get_gsn(p, x1, y) && before_igs == pict_get_gsn(p, x2, y) && before_igs == pict_get_gsn(p, x, y1) && before_igs == pict_get_gsn(p, x, y2) && before_igs == pict_get_gsn(p, x1, y1) && before_igs == pict_get_gsn(p, x2, y1) && before_igs == pict_get_gsn(p, x1, y2) && before_igs == pict_get_gsn(p, x2, y2) ) {
                                            checkpixels = 0;
                                             actual_igs = before_igs;
                                             }  }


					if (pict_is_validpixel(p, x, y) && before_igs > 0 && checkpixels==1) {
						actual_igs =
							find_near_pgs(p, x, y, max_angle, before_igs,
										  igs, 1);
						if (!(actual_igs == before_igs)) {
							change = 1;
						}
						if (before_igs > 0) {
							actual_igs = pict_get_gsn(p, x, y);
						/*	icol = setglcolor(p, x, y, actual_igs, uniform_gs, u_r, u_g , u_b); */
						}

					}
				}
			}
/*	pict_write(p,"secondscan.pic");*/
	}}

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
							                        pict_get_z1_gsn(p,igs) = 0;
							                        pict_get_z2_gsn(p,igs) = 0;

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
							                        pict_get_z1_gsn(p,igs) = 0;
							                        pict_get_z2_gsn(p,igs) = 0;

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
 /*											icol =
												setglcolor(p, x, y,
														   actual_igs, uniform_gs, u_r, u_g , u_b);*/
										}

									}
								}
							}

						}
				}


			}
		}
	}

/* calculation of direct vertical illuminance for CGI and for disability glare, coloring glare sources*/

	if (calcfast == 0) {
	for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
					if (pict_get_gsn(p, x, y) > 0) {
					        actual_igs = pict_get_gsn(p, x, y);
					        delta_E=DOT(p->view.vdir, pict_get_cached_dir(p, x, y))* pict_get_omega(p, x, y)* luminance(pict_get_color(p, x, y));
						pict_get_Eglare(p,actual_igs ) = pict_get_Eglare(p,actual_igs ) + delta_E;
						E_v_dir = E_v_dir + delta_E;
						setglcolor(p, x, y, actual_igs, uniform_gs, u_r, u_g , u_b);
					}
				}
			}
		}
	lum_backg_cos = (E_v - E_v_dir) / 3.1415927;

         }
/* calc of band luminance distribution if applied */
	if (band_calc == 1) {
		x_max = pict_get_xsize(p) - 1;
		y_max = pict_get_ysize(p) - 1;
		y_mid = (int)(y_max/2);
		for (x = 0; x <= x_max; x++)
		for (y = 0; y <= y_max; y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {

			r_actual = acos(DOT(pict_get_cached_dir(p, x, y_mid), pict_get_cached_dir(p, x, y))) ;
			act_lum = luminance(pict_get_color(p, x, y));

			if ((r_actual <= band_angle) || (y == y_mid) ) {
				lum_band[0]=luminance(pict_get_color(p, x, y));
	                        muc_rvar_add_sample(s_band, lum_band);
				act_lum = luminance(pict_get_color(p, x, y));
				lum_band_av += pict_get_omega(p, x, y) * act_lum;
				omega_band += pict_get_omega(p, x, y);
				if (band_color == 1) {
					pict_get_color(p, x, y)[RED] = 0.0;
					pict_get_color(p, x, y)[GRN] = act_lum / WHTEFFICACY / CIE_gf;
					pict_get_color(p, x, y)[BLU] = 0.0;
				}
			}
		}}}
                lum_band_av=lum_band_av/omega_band;
		muc_rvar_get_vx(s_band,lum_band_std);
	        muc_rvar_get_percentile(s_band,lum_band_median,0.75);
		per_75_band=lum_band_median[0];
	        muc_rvar_get_percentile(s_band,lum_band_median,0.95);
		per_95_band=lum_band_median[0];
	        muc_rvar_get_median(s_band,lum_band_median);
                muc_rvar_get_bounding_box(s_band,bbox_band);
	
 	        printf ("band:band_omega,band_av_lum,band_median_lum,band_std_lum,band_perc_75,band_perc_95,band_lum_min,band_lum_max: %f %f %f %f %f %f %f %f\n",omega_band,lum_band_av,lum_band_median[0],sqrt(lum_band_std[0]),per_75_band,per_95_band,bbox_band[0],bbox_band[1] );
/*	av_lum = av_lum / omega_sum; */
/*    printf("average luminance of band %f \n",av_lum);*/
/*      printf("solid angle of band %f \n",omega_sum);*/
	} 

/*printf("total number of glare sources: %i\n",igs); */
	lum_sources = 0;
	omega_sources = 0;
	i=0;
	for (x = 0; x <= igs; x++) {
	if (pict_get_npix(p, x) > 0) {
		lum_sources += pict_get_av_lum(p, x) * pict_get_av_omega(p, x);
		omega_sources += pict_get_av_omega(p, x);
		i=i+1;
	}}
       
        if (sang == omega_sources) {
               lum_backg =avlum;
        } else {
	       lum_backg =(sang * avlum - lum_sources) / (sang - omega_sources);
	}


        if (i == 0) {
               lum_sources =0.0;
        } else { lum_sources=lum_sources/omega_sources;
	}



	if (non_cos_lb == 0) {
			lum_backg = lum_backg_cos;
	}

/* file writing NOT here 
	if (checkfile == 1) {
		pict_write(p, file_out);
	}
*/

/* print detailed output */
	
	
	if (detail_out == 1) {

/* masking */

	if ( masking == 1 ) {
	
	
		for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {
                                        if (luminance(pict_get_color(pm, x, y))>0.1) {
/*					   printf ("hallo %i %i %f \n",x,y,luminance(pict_get_color(pm, x, y)));*/
					          i_mask=i_mask+1;
						  lum_mask[0]=luminance(pict_get_color(p, x, y));
						 /* printf ("%f \n",lum_mask[0]);*/
	                                          muc_rvar_add_sample(s_mask, lum_mask);
						  omega_mask +=	pict_get_omega(p, x, y);
						  lum_mask_av += pict_get_omega(p, x, y)* luminance(pict_get_color(p, x, y));
						  E_v_mask +=DOT(p->view.vdir, pict_get_cached_dir(p, x, y))*pict_get_omega(p, x, y)* luminance(pict_get_color(p, x, y));

						  pict_get_color(pm, x, y)[RED] = luminance(pict_get_color(p, x, y))/179.0;
					          pict_get_color(pm, x, y)[GRN] = luminance(pict_get_color(p, x, y))/179.0;
					          pict_get_color(pm, x, y)[BLU] = luminance(pict_get_color(p, x, y))/179.0;
  
                                           } else {
					   pict_get_color(p, x, y)[RED] = 0.0 ;
					   pict_get_color(p, x, y)[GRN] = 0.0 ;
					   pict_get_color(p, x, y)[BLU] = 0.0 ;
					   }
				}
			}
		}
/* Average luminance in masking area (weighted by solid angle) */
                lum_mask_av=lum_mask_av/omega_mask;
		muc_rvar_get_vx(s_mask,lum_mask_std);
	        muc_rvar_get_percentile(s_mask,lum_mask_median,0.75);
		per_75_mask=lum_mask_median[0];
	        muc_rvar_get_percentile(s_mask,lum_mask_median,0.95);
		per_95_mask=lum_mask_median[0];
	        muc_rvar_get_median(s_mask,lum_mask_median);
                muc_rvar_get_bounding_box(s_mask,bbox);
/* PSGV only why masking of window is applied! */
                 pgsv = get_pgsv(E_v, E_v_mask, omega_mask, lum_mask_av);
                 pgsv_sat =get_pgsv_sat(E_v);
	        printf ("masking:no_pixels,omega,av_lum,median_lum,std_lum,perc_75,perc_95,lum_min,lum_max,pgsv,pgsv_sat: %i %f %f %f %f %f %f %f %f %f %f\n",i_mask,omega_mask,lum_mask_av,lum_mask_median[0],sqrt(lum_mask_std[0]),per_75_mask,per_95_mask,bbox[0],bbox[1],pgsv,pgsv_sat );

		
		
	}

/* zones */

	if ( zones > 0 ) {
	
		for (x = 0; x < pict_get_xsize(p); x++)
		for (y = 0; y < pict_get_ysize(p); y++) {
			if (pict_get_hangle(p, x, y, p->view.vdir, p->view.vup, &ang)) {
				if (pict_is_validpixel(p, x, y)) {


			         r_actual =acos(DOT(pict_get_cached_dir(p, x, y), pict_get_cached_dir(p, x_zone, y_zone))) * 2;
			         lum_actual = luminance(pict_get_color(p, x, y));
			         act_gsn=pict_get_gsn(p, x, y);
			    /*     printf("1act_gsn,pict_get_z1_gsn,pict_get_z2_gsn %i %f %f \n",act_gsn,pict_get_z1_gsn(p,act_gsn),pict_get_z2_gsn(p,act_gsn));*/
			         
/*zone1 */
			if (r_actual <= angle_z1) {
					          i_z1=i_z1+1;
						  lum_z1[0]=lum_actual;
	                                          muc_rvar_add_sample(s_z1, lum_z1);
						  omega_z1 += pict_get_omega(p, x, y);
						  lum_z1_av += pict_get_omega(p, x, y)* lum_actual;
						  setglcolor(p,x,y,1,1 , 0.66, 0.01 ,0.33);
/*check if separation gsn already exist */

						   if (act_gsn > 0){
						   if (pict_get_z1_gsn(p,act_gsn) == 0) {
											  pict_new_gli(p);
											  igs = igs + 1;
							                        	  pict_get_z1_gsn(p,act_gsn) = igs*1.000;
							                        	  pict_get_z1_gsn(p,igs) = -1.0;
							                        	  pict_get_z2_gsn(p,igs) = -1.0;
							                        	  splitgs=(int)(pict_get_z1_gsn(p,act_gsn));
/*							                        	  printf ("Glare source %i gets glare source %i in zone 1 : %i %f %f \n",act_gsn,igs,splitgs,pict_get_z1_gsn(p,act_gsn),pict_get_z1_gsn(p,igs)); */ 
							                                 }
						    splitgs=(int)(pict_get_z1_gsn(p,act_gsn)); 
					/*	    printf("splitgs%i \n",splitgs);       */               
						    split_pixel_from_gs(p, x, y, splitgs, uniform_gs, u_r, u_g , u_b);
						}	                         
						  }
/*zone2 */

			if (r_actual > angle_z1 && r_actual<= angle_z2 ) {

					          i_z2=i_z2+1;
						  lum_z2[0]=lum_actual;
	                                          muc_rvar_add_sample(s_z2, lum_z2);
						  omega_z2 +=	pict_get_omega(p, x, y);
						  lum_z2_av += pict_get_omega(p, x, y)* lum_actual;
						  setglcolor(p,x,y,1,1 , 0.65, 0.33 ,0.02);
/*						  printf("zone 2 x y act_gsn pict_get_z1_gsn(p,act_gsn) pict_get_z2_gsn(p,act_gsn): %i %i %f 1:%f 2:%f \n",x,y,act_gsn,pict_get_z1_gsn(p,act_gsn),pict_get_z2_gsn(p,act_gsn));
*/						  if (act_gsn > 0){
						   if (pict_get_z2_gsn(p,act_gsn) == 0) {
											  pict_new_gli(p);
											  igs = igs + 1;
							                        	  pict_get_z2_gsn(p,act_gsn) = igs*1.000;
							                        	  pict_get_z1_gsn(p,igs) = -2.0;
							                        	  pict_get_z2_gsn(p,igs) = -2.0;
/*							                        	   printf ("Glare source %i gets glare source %i in zone 2 \n",act_gsn,igs ); */
							                                 }
						splitgs=(int)(pict_get_z2_gsn(p,act_gsn));
/*						printf("splitgs %i \n",splitgs);*/  	                        
						    split_pixel_from_gs(p, x, y, splitgs, uniform_gs, u_r, u_g , u_b);
                                           }
				}

			}
		} }
/* Average luminance in zones (weighted by solid angle), calculate percentiles, median min and max in zones  */
               if (zones == 2 ) {
                lum_z2_av=lum_z2_av/omega_z2;
		muc_rvar_get_vx(s_z2,lum_z2_std);
	        muc_rvar_get_percentile(s_z2,lum_z2_median,0.75);
		per_75_z2=lum_z2_median[0];
	        muc_rvar_get_percentile(s_z2,lum_z2_median,0.95);
		per_95_z2=lum_z2_median[0];
	        muc_rvar_get_median(s_z2,lum_z2_median);
                muc_rvar_get_bounding_box(s_z2,bbox_z2);
                }
                lum_z1_av=lum_z1_av/omega_z1;
 		muc_rvar_get_vx(s_z1,lum_z1_std);
	        muc_rvar_get_percentile(s_z1,lum_z1_median,0.75);
		per_75_z1=lum_z1_median[0];
	        muc_rvar_get_percentile(s_z1,lum_z1_median,0.95);
		per_95_z1=lum_z1_median[0];
	        muc_rvar_get_median(s_z1,lum_z1_median);
                muc_rvar_get_bounding_box(s_z1,bbox_z1);
 
 	        printf ("zoning:z1_omega,z1_av_lum,z1_median_lum,z1_std_lum,z1_perc_75,z1_perc_95,z1_lum_min,z1_lum_max: %f %f %f %f %f %f %f %f\n",omega_z1,lum_z1_av,lum_z1_median[0],sqrt(lum_z1_std[0]),per_75_z1,per_95_z1,bbox_z1[0],bbox_z1[1] );

               if (zones == 2 ) {

   	        printf ("zoning:z2_omega,z2_av_lum,z2_median_lum,z2_std_lum,z2_perc_75,z2_perc_95,z2_lum_min,z2_lum_max:  %f %f %f %f %f %f %f %f\n",omega_z2,lum_z2_av,lum_z2_median[0],sqrt(lum_z2_std[0]),per_75_z2,per_95_z2,bbox_z2[0],bbox_z2[1] );
 }             
		
	}

		i = 0;
		for (x = 0; x <= igs; x++) {
/* resorting glare source numbers */
			if (pict_get_npix(p, x) > 0) {
				i = i + 1;
			}
		}
		no_glaresources=i;

/* glare sources */
		printf
			("%i No pixels x-pos y-pos L_s Omega_s Posindx L_b L_t E_vert Edir Max_Lum Sigma xdir ydir zdir Eglare_cie Lveil_cie teta glare_zone\n",
			 i);
		if (i == 0) {
			printf("%i %f %f %f %f %.10f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", i, 0.0, 0.0,
				   0.0, 0.0, 0.0, 0.0, lum_backg, lum_task, E_v, E_v_dir,abs_max,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

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
					if (pict_get_z1_gsn(p,x)<0) {
					act_gsn=(int)(-pict_get_z1_gsn(p,x));
					}else{
					act_gsn=0;
					}
					printf("%i %f %f %f %f %.10f %f %f %f %f %f %f %f %f %f %f %f %f %f %i \n",
						   i, pict_get_npix(p, x), pict_get_av_posx(p, x),
						   pict_get_ysize(p) - pict_get_av_posy(p, x),
						   pict_get_av_lum(p, x), pict_get_av_omega(p, x),
						   get_posindex(p, pict_get_av_posx(p, x),
										pict_get_av_posy(p, x),
										posindex_2), lum_backg, lum_task,
						   E_v, E_v_dir, abs_max, sigma * 180 / 3.1415927
						   ,pict_get_cached_dir(p, x2,y2)[0],pict_get_cached_dir(p, x2,y2)[1],pict_get_cached_dir(p, x2,y2)[2],pict_get_Eglare(p,x),Lveil_cie,teta,act_gsn );
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
        lum_a= E_v2/3.1415927;
	dgi = get_dgi(p, lum_backg, igs, posindex_2);
	ugr = get_ugr(p, lum_backg, igs, posindex_2);
	ugp = get_ugp (p,lum_backg, igs, posindex_2);
	ugr_exp = get_ugr_exp (p,lum_backg_cos,lum_a, igs, posindex_2);
        dgi_mod = get_dgi_mod(p, lum_a, igs, posindex_2);
	cgi = get_cgi(p, E_v, E_v_dir, igs, posindex_2);
	dgr = get_dgr(p, avlum, igs, posindex_2);
	vcp = get_vcp(dgr);
	Lveil = get_disability(p, avlum, igs);
	if (no_glaresources == 0) {
 		dgi = 0.0;
		ugr = 0.0;
		ugp = 0.0;
		ugr_exp = 0.0;
		dgi_mod = 0.0;
		cgi = 0.0;
		dgr = 0.0;
		vcp = 100.0;
  	}
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
       	        muc_rvar_get_median(s_noposweight,lum_nopos_median);
       	        muc_rvar_get_median(s_posweight,lum_pos_median);
       	        muc_rvar_get_median(s_posweight2,lum_pos2_median);
 


        
			printf
				("dgp,av_lum,E_v,lum_backg,E_v_dir,dgi,ugr,vcp,cgi,lum_sources,omega_sources,Lveil,Lveil_cie,dgr,ugp,ugr_exp,dgi_mod,av_lum_pos,av_lum_pos2,med_lum,med_lum_pos,med_lum_pos2: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
				 dgp, avlum, E_v, lum_backg, E_v_dir, dgi, ugr, vcp, cgi,
				 lum_sources, omega_sources, Lveil,Lveil_cie_sum,dgr,ugp,ugr_exp,dgi_mod,lum_pos_mean,lum_pos2_mean/sang,lum_nopos_median[0],lum_pos_median[0],lum_pos2_median[0]);
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
