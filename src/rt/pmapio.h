/* RCSid $Id: pmapio.h,v 2.6 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ======================================================================
   Photon map file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapio.h,v 2.6 2016/05/17 17:39:47 rschregle Exp $
*/


#ifndef PMAPIO_H
   #define PMAPIO_H
   
   #include "pmapdata.h"


   /* File format version with feature suffix */
   #define PMAP_FILEVER_MAJ   "3.0"
   
#ifdef PMAP_OOC
   #define PMAP_FILEVER_TYP   "o"
#else
   #define PMAP_FILEVER_TYP   "k"
#endif

#ifdef PMAP_FLOAT_FLUX
   #define PMAP_FILEVER_FLX   "f"
#else
   #define PMAP_FILEVER_FLX   ""
#endif

#ifdef PMAP_PRIMARYPOS
   #define PMAP_FILEVER_PRI   "p"
#else
   #define PMAP_FILEVER_PRI   ""
#endif   
   
   #define PMAP_FILEVER       (PMAP_FILEVER_MAJ PMAP_FILEVER_TYP \
                               PMAP_FILEVER_FLX PMAP_FILEVER_PRI)
   
   
   void savePhotonMap (const PhotonMap *pmap, const char *fname,
                       int argc, char **argv);
   /* Save portable photon map of specified type to filename,
      along with the corresponding command line. */

   PhotonMapType loadPhotonMap (PhotonMap *pmap, const char* fname);
   /* Load portable photon map from specified filename, and return photon
    * map type as identified from file header. */

#endif
