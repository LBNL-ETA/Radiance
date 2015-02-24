/* 
   ==================================================================
   Photon map interface to RADIANCE raycalls

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================   
   
   $Id: pmapray.h,v 2.1 2015/02/24 19:39:27 greg Exp $
*/


#include "ray.h"


void ray_init_pmap ();
/* Interface to ray_init() and rtmain/rpmain/rvmain; init & load pmaps */

void ray_done_pmap ();
/* Interface to ray_done() and rtmain/rpmain/rvmain; free photon maps */

void ray_save_pmap (RAYPARAMS *rp);
/* Interface to ray_save(); save photon map params */

void ray_restore_pmap (RAYPARAMS *rp);
/* Interface to ray_restore(); restore photon mapping params */
