/* 
   ==================================================================
   Photon map interface to RADIANCE ambient calculation
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: pmapamb.h,v 2.1 2015/02/24 19:39:26 greg Exp $
*/


#ifndef PMAPAMB_H
#define PMAPAMB_H

   #include "pmapdata.h"   

   int ambPmap (COLOR aval, RAY *r, int rdepth);
   /* Factor irradiance from global photon map into ambient coefficient
    * aval; return 1 on success, else 0 (with aval unmodified) */

   int ambPmapCaustic (COLOR aval, RAY *r, int rdepth);
   /* Factor irradiance from caustic photon map into ambient coeffiecient
    * aval; return 1 if successful, else 0 (with aval set to zero) */
#endif
