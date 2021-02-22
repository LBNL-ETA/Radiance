/* RCSid $Id: pmapmat.h,v 2.14 2021/02/22 13:27:49 rschregle Exp $ */
/* 
   ======================================================================
   Photon map support routines for scattering by materials. 
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme FARESYS")
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ======================================================================
   
*/


#ifndef PMAPMAT_H
   #define PMAPMAT_H

   #include "pmap.h"

   /* 
      Check for paths already accounted for in photon map to avoid
      double-counting during backward raytracing.
      
      ambRayInPmap():      Check for DIFFUSE -> (DIFFUSE|SPECULAR) -> *	 
                           subpaths.  These occur during the backward pass	 
                           when an ambient ray spawns an indirect diffuse or	 
                           specular ray.  These subpaths are already	 
                           accounted for in the global photon map.
   */
   #define ambRayInPmap(r) ((r) -> crtype & AMBIENT && photonMapping && \
                              (ambounce < 0 || (r) -> rlvl > 1))

   /* 
      shadowRayInPmap():   Check for DIFFUSE -> SPECULAR -> LIGHT
                           subpaths. These occur during the backward pass 
                           when a shadow ray is transferred through a
                           transparent material. These subpaths are already
                           accounted for by caustic photons in the global,
                           contrib, or dedicated caustic photon map.
      
      !!! DISABLED FOR TESTING PENDING REPLACEMENT BY srcRayInPmap() !!!                           
   */
#if 1
   #define shadowRayInPmap(r) 0
#else   
   #define shadowRayInPmap(r) (((globalPmap && ambounce < 0) || \
                                 causticPmap || contribPmap) && \
                                 (r) -> crtype & SHADOW)
#endif

   /* srcRayInPmap():      Check whether a source ray transferred through
    *                      medium (e.g.  glass/dielectric) is already
    *                      accounted for in the photon map.  This is used by
    *                      badcomponent() in source.c when checking source
    *                      hits (glows and lights, hence ambient and shadow
    *                      rays).
    */
   #define srcRayInPmap(r)    (((globalPmap && ambounce < 0) || \
                                 causticPmap || contribPmap) && \
                                 (r) -> crtype & (AMBIENT | SHADOW) && \
                                 (r) -> rtype & (TRANS | REFRACTED))

   /* Check if scattered ray spawns a caustic photon; 
    * !!! NOTE this returns a single bit as boolean value (0|1), rather
    * !!! than the original short int, hence the explicit test against zero.
    * !!! This allows the macro the be used in a conditional statement
    * !!! and when setting a photon's caustic flag in newPhoton(). */
   #define PMAP_CAUSTICRAY(r)    (((r) -> rtype & SPECULAR) != 0) 


   /* Scattered photon ray types for photonRay() */
   #define  PMAP_DIFFREFL        (REFLECTED | AMBIENT)
   #define  PMAP_DIFFTRANS       (REFLECTED | AMBIENT | TRANS)
   #define  PMAP_SPECREFL        (REFLECTED | SPECULAR)
   #define  PMAP_SPECTRANS       (REFLECTED | SPECULAR | TRANS)
   #define  PMAP_REFRACT         (REFRACTED | SPECULAR)
   #define  PMAP_XFER            (TRANS)



   /* Dispatch table for photon scattering functions */
   extern int (*photonScatter []) (OBJREC*, RAY*);
   
   /* List of antimatter sensor modifier names */
   extern char *photonSensorList [MAXSET + 1];



   /* Spawn a new photon ray from a previous one; this is effectively a
    * customised rayorigin(). */
   void photonRay (const RAY *rayIn, RAY *rayOut, int rayOutType, 
                   COLOR fluxAtten);

   /* Init photonScatter[] dispatch table with material specific scattering
      routines */
   void initPhotonScatterFuncs ();

   /* Find antimatter geometry declared as photon sensors */   
   void getPhotonSensors (char **sensorList);
   
#endif
