/* RCSid $Id: pmapray.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ==================================================================
   Photon map interface to RADIANCE raycalls

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapray.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $
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
