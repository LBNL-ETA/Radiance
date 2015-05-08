/* 
   ==================================================================
   Photon map file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapio.h,v 2.2 2015/05/08 13:20:22 rschregle Exp $
*/


#ifndef PMAPIO_H
   #define PMAPIO_H
   
   #include "pmapdata.h"
   
   #define PMAP_FILEVER 160u              /* File format version */
   
   void savePhotonMap (const PhotonMap *pmap, const char *fname,
                       PhotonMapType type, int argc, char **argv);
   /* Save portable photon map of specified type to filename,
      along with the corresponding command line. */

   PhotonMapType loadPhotonMap (PhotonMap *pmap, const char* fname);
   /* Load portable photon map from specified filename, and return photon
    * map type as identified from file header. */

#endif
