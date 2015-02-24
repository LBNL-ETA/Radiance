/* 
   ==================================================================
   Photon map types and corresponding file format strings

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: pmaptype.c,v 2.1 2015/02/24 19:39:27 greg Exp $
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
