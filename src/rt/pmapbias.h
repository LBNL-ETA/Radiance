/* RCSid $Id: pmapbias.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ==================================================================
   Bias compensation for photon density estimates
   
   For background see:
      R. Schregle, "Bias Compensation for Photon Maps",
      Computer Graphics Forum, v22:n4, pp. 729-742, Dec. 2003.

   Roland Schregle (roland.schregle@gmail.com)
   (c) Fraunhofer Institute for Solar Energy Systems
   ==================================================================
      
   $Id: pmapbias.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $
*/


#ifndef PMAPBIASCOMP_H
   #define PMAPBIASCOMP_H

   #include "pmapdata.h"
      
   /* Bias compensation weighting function */
   /* #define BIASCOMP_WGT(n) 1 */
   /* #define BIASCOMP_WGT(n) (n) */
   #define BIASCOMP_WGT(n) ((n) * (n))
   /* #define BIASCOMP_WGT(n) ((n) * (n) * (n)) */
   /* #define BIASCOMP_WGT(n) exp(0.003 * (n)) */

   /* Dump photon bandwidth for bias compensated density estimates */
   /* #define BIASCOMP_BWIDTH */

   void biasComp (PhotonMap*, COLOR);
   /* Photon density estimate with bias compensation, returning irradiance. 
      Expects photons in search queue after a kd-tree lookup. */

   void volumeBiasComp (PhotonMap*, const RAY*, COLOR);   
   /* Photon volume density estimate with bias compensation, returning
      irradiance. Expects photons in search queue after a kd-tree lookup. */

#endif
