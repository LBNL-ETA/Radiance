#ifndef lint
static const char RCSid[] = "$Id: pmapbias.c,v 2.5 2016/05/17 17:39:47 rschregle Exp $";
#endif

/* 
   ==================================================================
   Bias compensation for photon density estimates
   
   For background see:
      R. Schregle, "Bias Compensation for Photon Maps",
      Computer Graphics Forum, v22:n4, pp. 729-742, Dec. 2003.

   Roland Schregle (roland.schregle@gmail.com)
   (c) Fraunhofer Institute for Solar Energy Systems
   ==================================================================
   
   $Id: pmapbias.c,v 2.5 2016/05/17 17:39:47 rschregle Exp $    
*/



#include "pmapbias.h"
#include "pmap.h"
#include "pmaprand.h"



void squeuePartition (PhotonSearchQueueNode* squeue, unsigned lo, 
                      unsigned mid, unsigned hi)
/* REVERSE Partition squeue such that all photons in 
   squeue-hi+1..squeue-mid are farther than the median at squeue-mid+1, 
   and those in squeue-mid+2..squeue-lo+1 are closer than the median. 
   This means that squeue points to the END of the queue, and the (1-based) 
   indices are offsets relative to it. This convoluted scheme is adopted 
   since the queue is initially a maxheap, so reverse sorting is expected 
   to be faster. */
{
   unsigned                l, h, p;
   PhotonSearchQueueNode   *lp, *hp, *pp, tmp;
   float                   pivot;

   while (hi > lo) {
      /* Grab pivot node in middle as an educated guess, since our
         queue is sorta sorted. */
      l = lo;
      h = hi;
      p = mid;
      lp = squeue - lo + 1;
      hp = squeue - hi + 1;
      pp = squeue - p + 1;
      pivot = pp -> dist2;
      
      /* l & h converge, swapping elements out of order with respect to 
         pivot node. */
      while (l < h) {
         while (lp -> dist2 <= pivot && l <= h && l < hi)
            ++l, --lp;
            
         while (hp -> dist2 >= pivot && h >= l && h > lo)
            --h, ++hp;
            
         if (l < h) {
            /* Swap */
            tmp.idx = lp -> idx;
            tmp.dist2 = lp -> dist2;
            lp -> idx = hp -> idx;
            lp -> dist2 = hp -> dist2;
            hp -> idx = tmp.idx;
            hp -> dist2 = tmp.dist2;
         }
      }
      
      /* Swap convergence and pivot node */
      if (p > h) {
         /* Need this otherwise shit happens! 
            Since lp -> dist2 > hp -> dist2, we swap either l or p depending 
            on whether we're above or below p */
         h = l;
         hp = lp;
      }
      
      tmp.idx = hp -> idx;
      tmp.dist2 = hp -> dist2;
      hp -> idx = pp -> idx;
      hp -> dist2 = pivot;
      pp -> idx = tmp.idx;
      pp -> dist2 = tmp.dist2; 
      
      if (h >= mid) 
         hi = h - 1; 
         
      if (h <= mid) 
         lo = h + 1; 
   }
   
   /* Once lo & hi have converged, we have found the median! */
}



void biasComp (PhotonMap* pmap, COLOR irrad)
/* Photon density estimate with bias compensation -- czech dis shit out! */
{
   unsigned                i, numLo, numHi, numMid;
   float                   r, totalWeight = 0;
   COLOR                   fluxLo, fluxMid, irradVar, irradAvg, p, d;
   Photon                  *photon;
   PhotonSearchQueueNode   *sqn, *sqEnd;
   PhotonBiasCompNode      *hist, *histEnd;
   
   if (!pmap -> biasCompHist) {
      /* Allocate bias compensation history */
      numHi = pmap -> maxGather - pmap -> minGather;
      for (i = pmap -> minGather + 1; numHi > 1; numHi >>= 1, ++i);      
      pmap -> biasCompHist = calloc(i, sizeof(PhotonBiasCompNode));      
       
      if (!pmap -> biasCompHist)
         error(USER, "can't allocate bias compensation history");
   }
   
   numLo = min(pmap -> minGather, pmap -> squeue.tail - 1);
   numHi = min(pmap -> maxGather, pmap -> squeue.tail - 1);
   sqn = sqEnd = pmap -> squeue.node + pmap -> squeue.tail - 1;
   histEnd = pmap -> biasCompHist;
   setcolor(fluxLo, 0, 0, 0);
   
   /* Partition to get numLo closest photons starting from END of queue */
   squeuePartition(sqEnd, 1, numLo + 1, numHi);
   
   /* Get irradiance estimates (ignoring 1 / PI) using 1..numLo photons 
      and chuck in history. Queue is traversed BACKWARDS. */
   for (i = 1; i <= numLo; i++, sqn--) {
      /* Make sure furthest two photons are consecutive w.r.t. distance */
      squeuePartition(sqEnd, i, i + 1, numLo + 1);
      photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
      getPhotonFlux(photon, irrad);
      addcolor(fluxLo, irrad);
      /* Average radius between furthest two photons to improve accuracy */
      r = 0.25 * (sqn -> dist2 + (sqn - 1) -> dist2 +
                  2 * sqrt(sqn -> dist2 * (sqn - 1) -> dist2));
      /* Add irradiance and weight to history. Weights should increase
         monotonically based on number of photons used for the estimate. */
      histEnd -> irrad [0] = fluxLo [0] / r;
      histEnd -> irrad [1] = fluxLo [1] / r;
      histEnd -> irrad [2] = fluxLo [2] / r;
      totalWeight += histEnd++ -> weight = BIASCOMP_WGT((float)i);
   }
   
   /* Compute expected value (average) and variance of irradiance based on 
      history for numLo photons. */
   setcolor(irradAvg, 0, 0, 0);
   setcolor(irradVar, 0, 0, 0);
   
   for (hist = pmap -> biasCompHist; hist < histEnd; ++hist)
      for (i = 0; i <= 2; ++i) {
         irradAvg [i] += r = hist -> weight * hist -> irrad [i];
         irradVar [i] += r * hist -> irrad [i];
      }
      
   for (i = 0; i <= 2; ++i) {
      r = irradAvg [i] /= totalWeight;
      irradVar [i] = irradVar [i] / totalWeight - r * r;
   }
   
   /* Do binary search within interval [numLo, numHi]. numLo is towards
      the END of the queue. */
   while (numHi - numLo > 1) {
      numMid = (numLo + numHi) >> 1;
      /* Split interval to get numMid closest photons starting from the 
         END of the queue */
      squeuePartition(sqEnd, numLo, numMid, numHi);
      /* Make sure furthest two photons are consecutive wrt distance */
      squeuePartition(sqEnd, numMid, numMid + 1, numHi);
      copycolor(fluxMid, fluxLo);
      sqn = sqEnd - numLo;
      
      /* Get irradiance for numMid photons (ignoring 1 / PI) */
      for (i = numLo; i < numMid; i++, sqn--) {
         photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
         getPhotonFlux(photon, irrad);
         addcolor(fluxMid, irrad);
      }
      
      /* Average radius between furthest two photons to improve accuracy */
      r = 0.25 * (sqn -> dist2 + (sqn + 1) -> dist2 +
                  2 * sqrt(sqn -> dist2 * (sqn + 1) -> dist2));
                  
      /* Get deviation from current average, and probability that it's due
         to noise from gaussian distribution based on current variance. Since
         we are doing this for each colour channel we can also detect
         chromatic bias. */
      for (i = 0; i <= 2; ++i) {
         d [i] = irradAvg [i] - (irrad [i] = fluxMid [i] / r);
         p [i] = exp(-0.5 * d [i] * d [i] / irradVar [i]);
      }
      
      if (pmapRandom(pmap -> randState) < colorAvg(p)) {
         /* Deviation is probably noise, so add mid irradiance to history */
         copycolor(histEnd -> irrad, irrad);
         totalWeight += histEnd++ -> weight = BIASCOMP_WGT((float)numMid);
         setcolor(irradAvg, 0, 0, 0);
         setcolor(irradVar, 0, 0, 0);
         
         /* Update average and variance */
         for (hist = pmap -> biasCompHist; hist < histEnd; ++hist)
            for (i = 0; i <= 2; i++) {
               r = hist -> irrad [i];
               irradAvg [i] += hist -> weight * r;
               irradVar [i] += hist -> weight * r * r;
            }
            
         for (i = 0; i <= 2; i++) {
            r = irradAvg [i] /= totalWeight;
            irradVar [i] = irradVar [i] / totalWeight - r * r;
         }
         
         /* Need more photons --> recurse on [numMid, numHi] */
         numLo = numMid;
         copycolor(fluxLo, fluxMid);
      }
      else 
         /* Deviation is probably bias --> need fewer photons, 
            so recurse on [numLo, numMid] */
         numHi = numMid;
   }   

   --histEnd;
   
   for (i = 0; i <= 2; i++) {
      /* Estimated relative error */
      d [i] = histEnd -> irrad [i] / irradAvg [i] - 1;
      
#ifdef BIASCOMP_BWIDTH
      /* Return bandwidth instead of irradiance */
      irrad [i] = numHi / PI;
#else
      /* 1 / PI required as ambient normalisation factor */
      irrad [i] = histEnd -> irrad [i] / (PI * PI);
#endif
   }
   
   /* Update statistix */
   r = colorAvg(d);
   pmap -> rmsError += r * r;
   
   if (r < pmap -> minError) 
      pmap -> minError = r;
      
   if (r > pmap -> maxError) 
      pmap -> maxError = r;
   
   if (numHi < pmap -> minGathered) 
      pmap -> minGathered = numHi;
            
   if (numHi > pmap -> maxGathered) 
      pmap -> maxGathered = numHi;
      
   pmap -> totalGathered += numHi;
   ++pmap -> numDensity;
}



/* Volume bias compensation disabled (probably redundant) */
#if 0
void volumeBiasComp (PhotonMap* pmap, const RAY* ray, COLOR irrad)
/* Photon volume density estimate with bias compensation -- czech dis 
   shit out! */
{
   unsigned i, numLo, numHi, numMid = 0;
   PhotonSQNode *sq;
   PhotonBCNode *hist;
   const float gecc2 = ray -> gecc * ray -> gecc;
   float r, totalWeight = 0;
   PhotonSQNode *squeueEnd;
   PhotonBCNode *histEnd;
   COLOR fluxLo, fluxMid, irradVar, irradAvg, p, d;

   if (!pmap -> biasCompHist) {
      /* Allocate bias compensation history */
      numHi = pmap -> maxGather - pmap -> minGather; 
      for (i = pmap -> minGather + 1; numHi > 1; numHi >>= 1, ++i);
      pmap -> biasCompHist = (PhotonBCNode*)malloc(i * sizeof(PhotonBCNode));
      if (!pmap -> biasCompHist)
         error(USER, "can't allocate bias compensation history");
   }
   
   numLo = min(pmap -> minGather, pmap -> squeueEnd - 1);
   numHi = min(pmap -> maxGather, pmap -> squeueEnd - 1);
   sq = squeueEnd = pmap -> squeue + pmap -> squeueEnd - 1;
   histEnd = pmap -> biasCompHist;
   setcolor(fluxLo, 0, 0, 0);
   /* Partition to get numLo closest photons starting from END of queue */
   squeuePartition(squeueEnd, 1, numLo, numHi);
   
   /* Get irradiance estimates (ignoring constants) using 1..numLo photons 
      and chuck in history. Queue is traversed BACKWARDS. */
   for (i = 1; i <= numLo; i++, sq--) {
      /* Make sure furthest two photons are consecutive wrt distance */
      squeuePartition(squeueEnd, i, i + 1, numHi);
      
      /* Compute phase function for inscattering from photon */
      if (gecc2 <= FTINY) 
         r = 1;
      else {
         r = DOT(ray -> rdir, sq -> photon -> norm) / 127;
         r = 1 + gecc2 - 2 * ray -> gecc * r;
         r = (1 - gecc2) / (r * sqrt(r));
      }
      
      getPhotonFlux(sq -> photon, irrad);
      scalecolor(irrad, r);
      addcolor(fluxLo, irrad);
      /* Average radius between furthest two photons to improve accuracy */
      r = 0.25 * (sq -> dist + (sq - 1) -> dist +
                  2 * sqrt(sq -> dist * (sq - 1) -> dist));  
      r *= sqrt(r);
      /* Add irradiance and weight to history. Weights should increase
         monotonically based on number of photons used for the estimate. */     
      histEnd -> irrad [0] = fluxLo [0] / r;
      histEnd -> irrad [1] = fluxLo [1] / r;
      histEnd -> irrad [2] = fluxLo [2] / r;
      totalWeight += histEnd++ -> weight = BIASCOMP_WGT((float)i);
   }
   
   /* Compute expected value (average) and variance of irradiance based on 
      history for numLo photons. */
   setcolor(irradAvg, 0, 0, 0);
   setcolor(irradVar, 0, 0, 0);
   
   for (hist = pmap -> biasCompHist; hist < histEnd; ++hist)
      for (i = 0; i <= 2; ++i) {
         irradAvg [i] += r = hist -> weight * hist -> irrad [i];
         irradVar [i] += r * hist -> irrad [i];
      }
      
   for (i = 0; i <= 2; ++i) {
      r = irradAvg [i] /= totalWeight;
      irradVar [i] = irradVar [i] / totalWeight - r * r;
   }
   
   /* Do binary search within interval [numLo, numHi]. numLo is towards
      the END of the queue. */
   while (numHi - numLo > 1) {
      numMid = (numLo + numHi) >> 1;
      /* Split interval to get numMid closest photons starting from the 
         END of the queue */
      squeuePartition(squeueEnd, numLo, numMid, numHi);
      /* Make sure furthest two photons are consecutive wrt distance */
      squeuePartition(squeueEnd, numMid, numMid + 1, numHi);
      copycolor(fluxMid, fluxLo);
      sq = squeueEnd - numLo;
      
      /* Get irradiance for numMid photons (ignoring constants) */
      for (i = numLo; i < numMid; i++, sq--) {
         /* Compute phase function for inscattering from photon */
         if (gecc2 <= FTINY) 
            r = 1;
         else {
            r = DOT(ray -> rdir, sq -> photon -> norm) / 127;
            r = 1 + gecc2 - 2 * ray -> gecc * r;
            r = (1 - gecc2) / (r * sqrt(r));
         }
         
         getPhotonFlux(sq -> photon, irrad);
         scalecolor(irrad, r);
         addcolor(fluxMid, irrad);
      }
      
      /* Average radius between furthest two photons to improve accuracy */
      r = 0.25 * (sq -> dist + (sq + 1) -> dist +
                  2 * sqrt(sq -> dist * (sq + 1) -> dist));
      r *= sqrt(r);
      
      /* Get deviation from current average, and probability that it's due
         to noise from gaussian distribution based on current variance. Since
         we are doing this for each colour channel we can also detect
         chromatic bias. */
      for (i = 0; i <= 2; ++i) {
         d [i] = irradAvg [i] - (irrad [i] = fluxMid [i] / r);
         p [i] = exp(-0.5 * d [i] * d [i] / irradVar [i]);
      }
      
      if (pmapRandom(pmap -> randState) < colorAvg(p)) {
         /* Deviation is probably noise, so add mid irradiance to history */
         copycolor(histEnd -> irrad, irrad);
         totalWeight += histEnd++ -> weight = BIASCOMP_WGT((float)numMid);
         setcolor(irradAvg, 0, 0, 0);
         setcolor(irradVar, 0, 0, 0);
         
         /* Update average and variance */
         for (hist = pmap -> biasCompHist; hist < histEnd; ++hist)
            for (i = 0; i <= 2; i++) {
               r = hist -> irrad [i];
               irradAvg [i] += hist -> weight * r;
               irradVar [i] += hist -> weight * r * r;
            }
         for (i = 0; i <= 2; i++) {
            r = irradAvg [i] /= totalWeight;
            irradVar [i] = irradVar [i] / totalWeight - r * r;
         }
         
         /* Need more photons -- recurse on [numMid, numHi] */
         numLo = numMid;
         copycolor(fluxLo, fluxMid);
      }
      else 
         /* Deviation is probably bias -- need fewer photons, 
            so recurse on [numLo, numMid] */
         numHi = numMid;
   }   
   
   --histEnd;
   for (i = 0; i <= 2; i++) {
      /* Estimated relative error */
      d [i] = histEnd -> irrad [i] / irradAvg [i] - 1;
      /* Divide by 4 / 3 * PI for search volume (r^3 already accounted
         for) and phase function normalization factor 1 / (4 * PI) */
      irrad [i] = histEnd -> irrad [i] * 3 / (16 * PI * PI);
   }
   
   /* Update statistix */
   r = colorAvg(d);
   if (r < pmap -> minError) 
      pmap -> minError = r;
   if (r > pmap -> maxError) 
      pmap -> maxError = r;
   pmap -> rmsError += r * r;
   
   if (numMid < pmap -> minGathered) 
      pmap -> minGathered = numMid;
   if (numMid > pmap -> maxGathered) 
      pmap -> maxGathered = numMid;
      
   pmap -> totalGathered += numMid;
   ++pmap -> numDensity;
}
#endif
