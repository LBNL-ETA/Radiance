#ifndef lint
static const char RCSid[] = "$Id: pmaptype.c,v 2.5 2016/05/17 17:39:47 rschregle Exp $";
#endif

/* 
   ==================================================================
   Photon map types and corresponding file format strings

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmaptype.c,v 2.5 2016/05/17 17:39:47 rschregle Exp $
*/


#include "pmaptype.h"

#ifdef PMAP_OOC
   #define PMAP_FMTSUFFIX  "OOC_Photon_Map"
#else
   #define PMAP_FMTSUFFIX  "kdT_Photon_Map"
#endif


/* Format strings for photon map files corresponding to PhotonMapType */
const char *pmapFormat [NUM_PMAP_TYPES] = { 
   "Radiance_Global_"   PMAP_FMTSUFFIX, "Radiance_PreComp_" PMAP_FMTSUFFIX,
   "Radiance_Caustic_"  PMAP_FMTSUFFIX, "Radiance_Volume_"  PMAP_FMTSUFFIX,
   "Radiance_Direct_"   PMAP_FMTSUFFIX, "Radiance_Contrib_" PMAP_FMTSUFFIX
};

/* Photon map names per type */
const char *pmapName [NUM_PMAP_TYPES] = {
   "global", "precomp", "caustic", "volume", "direct", "contrib"
};
