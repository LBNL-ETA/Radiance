/* RCSid $Id: pmaptype.h,v 2.6 2024/09/17 16:36:05 greg Exp $ */

/* 
   ==================================================================
   Photon map types and corresponding file format strings

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmaptype.h,v 2.6 2024/09/17 16:36:05 greg Exp $
*/


#ifndef PMAPTYPE_H
   #define PMAPTYPE_H
   
#ifdef __cplusplus
extern "C" {
#endif

   /* Photon map types */
   typedef enum {
      PMAP_TYPE_NONE = -1, PMAP_TYPE_GLOBAL,    PMAP_TYPE_PRECOMP, 
      PMAP_TYPE_CAUSTIC,   PMAP_TYPE_VOLUME,    PMAP_TYPE_DIRECT,
      PMAP_TYPE_CONTRIB,   NUM_PMAP_TYPES
   } PhotonMapType;

   /* Check for valid photon map type */
   #define validPmapType(t)   ((t) >= 0 && (t) < NUM_PMAP_TYPES)
   
   /* Glob string for extracting photon map format from file header */
   #define PMAP_FORMAT_GLOB   "Radiance_*_Photon_Map"
   
   /* Format strings for photon map files corresponding to PhotonMapType */
   extern const char *pmapFormat [];      
   
   /* Photon map names per type */
   extern const char *pmapName [];

#ifdef __cplusplus
}
#endif

#endif
