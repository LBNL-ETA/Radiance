#ifndef lint
static const char RCSid[] = "$Id: pmaptype.c,v 2.3 2015/08/18 18:45:55 greg Exp $";
#endif
/* 
   ==================================================================
   Photon map types and corresponding file format strings

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmaptype.c,v 2.3 2015/08/18 18:45:55 greg Exp $
*/


#include "pmaptype.h"


/* Format strings for photon map files corresponding to PhotonMapType */
const char *pmapFormat [NUM_PMAP_TYPES] = { 
   "Radiance_Global_Photon_Map",    "Radiance_PreComp_Photon_Map",
   "Radiance_Caustic_Photon_Map",   "Radiance_Volume_Photon_Map",
   "Radiance_Direct_Photon_Map",    "Radiance_Contrib_Photon_Map"
};

/* Photon map names per type */
const char *pmapName [NUM_PMAP_TYPES] = {
   "global", "precomp", "caustic", "volume", "direct", "contrib"
};
