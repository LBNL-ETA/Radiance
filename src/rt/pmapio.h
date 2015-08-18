/* RCSid $Id: pmapio.h,v 2.4 2015/08/18 18:45:55 greg Exp $ */
/* 
   ==================================================================
   Photon map file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapio.h,v 2.4 2015/08/18 18:45:55 greg Exp $
*/


#ifndef PMAPIO_H
   #define PMAPIO_H
   
   #include "pmapdata.h"
   
   #define PMAP_FILEVER 160u              /* File format version */
   
   void savePhotonMap (const PhotonMap *pmap, const char *fname,
                       int argc, char **argv);
   /* Save portable photon map of specified type to filename,
      along with the corresponding command line. */

   PhotonMapType loadPhotonMap (PhotonMap *pmap, const char* fname);
   /* Load portable photon map from specified filename, and return photon
    * map type as identified from file header. */

#endif
