/* RCSid $Id: pmapray.h,v 2.8 2024/09/17 16:36:05 greg Exp $ */

/* 
   ==================================================================
   Photon map interface to RADIANCE raycalls

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapray.h,v 2.8 2024/09/17 16:36:05 greg Exp $
*/

/* Include after ray.h */

#ifdef __cplusplus
extern "C" {
#endif

void ray_init_pmap (void);
/* Interface to ray_init() and rtmain/rpmain/rvmain; init & load pmaps */

void ray_done_pmap (void);
/* Interface to ray_done() and rtmain/rpmain/rvmain; free photon maps */

void ray_save_pmap (RAYPARAMS *rp);
/* Interface to ray_save(); save photon map params */

void ray_restore_pmap (RAYPARAMS *rp);
/* Interface to ray_restore(); restore photon mapping params */

void ray_defaults_pmap (RAYPARAMS *rp);
/* Interface to ray_defaults(); set photon mapping defaults */

#ifdef __cplusplus
}
#endif
