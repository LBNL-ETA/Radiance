/* 
   ==================================================================
   Photon map file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: pmapio.h,v 2.1 2015/02/24 19:39:27 greg Exp $
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
