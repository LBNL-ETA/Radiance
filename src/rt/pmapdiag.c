#ifndef lint
static const char RCSid[] = "$Id: pmapdiag.c,v 2.8 2021/03/26 23:47:03 rschregle Exp $";
#endif

/* 
   ==================================================================
   Photon map diagnostic output and progress reports

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapdiag.c,v 2.8 2021/03/26 23:47:03 rschregle Exp $
*/



#include "pmap.h"
#include "pmapdiag.h"
#include "pmapdata.h"
#include "standard.h"



time_t repStartTime, repLastTime = 0;   /* Time at start & last report */
unsigned long repProgress,              /* Report progress counter */
              repComplete,              /* Report completion count */
              repEmitted;               /* Num emitted photons */


static char* biasCompStats (const PhotonMap *pmap, PhotonMapType type, 
                            char *stats)
/* Dump bias compensation statistics */
{
   unsigned avgBwidth;
   float avgErr;
   
   stats [0] = '\0';
   
   /* Check if photon map is valid and bias compensated */
   if (pmap && pmap -> maxGather > pmap -> minGather) {      
      avgBwidth = pmap -> numDensity 
                  ? pmap -> totalGathered / pmap -> numDensity : 0,
      avgErr = pmap -> numDensity 
               ? sqrt(pmap -> rmsError / pmap -> numDensity) : 0;
                                              
      sprintf(stats, "%d/%d/%d %s photon bwidth (%.1f/%.1f/%.1f%% error), ",
              pmap -> minGathered, pmap -> maxGathered, avgBwidth, 
              pmapName [type],
              100 * pmap -> minError, 100 * pmap -> maxError, 100 * avgErr);
                 
      return stats;
   }
   
   return NULL;
}



void pmapBiasCompReport (char *stats)
/* Append full bias compensation statistics to stats; interface to rpict's
 * report() */
{
   char tmp [512];
   unsigned t;
   
   stats [0] = '\0';
   
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (biasCompStats(photonMaps [t], t, tmp))
         strcat(stats, tmp);
}



#ifdef PMAP_OOC
   static char* oocCacheStats (const PhotonMap *pmap, PhotonMapType type, 
                               char *stats)
   /* Dump OOC I/O cache statistics */
   {
      const OOC_Cache   *cache;
      stats [0] = '\0';

      
      /* Check for photon map is valid and caching enabled */
      if (pmap && (cache = pmap -> store.cache) && cache -> numReads) {
         sprintf(stats, "%ld %s photons cached in %d/%d pages "
                 "(%.1f%% hit, %.1f%% rept, %.1f coll), ",
                 (unsigned long)cache -> pageCnt * cache -> recPerPage,
                 pmapName [type], cache -> pageCnt, cache -> numPages, 
                 100.0 * cache -> numHits / cache -> numReads,
                 100.0 * cache -> numRept / cache -> numReads,
                 (float)cache -> numColl / cache -> numReads);
                    
         return stats;
      }
      
      return NULL;
   }



   void pmapOOCCacheReport (char *stats)
   /* Append full OOC I/O cache statistics to stats; interface to rpict's
    * report() */
   {
      char tmp [512];
      unsigned t;
      
      stats [0] = '\0';
      
      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (oocCacheStats(photonMaps [t], t, tmp))
            strcat(stats, tmp);
   }
#endif



#ifndef NON_POSIX
   void pmapPreCompReport()
   /* Report global photon precomputation progress */
   {
      extern char *myhostname();
      float u, s, rtime, eta, progress;
      char tmp [512];
      
      #ifdef BSD
         struct rusage rubuf;
      #else
         struct tms tbuf;
         float period;
      #endif

      repLastTime = time(NULL);
      
      #ifdef BSD
         getrusage(RUSAGE_SELF, &rubuf);
         u = rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s = rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
         getrusage(RUSAGE_CHILDREN, &rubuf);
         u += rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s += rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
      #else
         times(&tbuf);

         #ifdef _SC_CLK_TCK
            period = 1.0 / sysconf(_SC_CLK_TCK);
         #else
            period = 1.0 / 60;
         #endif
      
         u = (tbuf.tms_utime + tbuf.tms_cutime) * period;
         s = (tbuf.tms_stime + tbuf.tms_cstime) * period;
      #endif

      sprintf(errmsg, "%lu precomputed, ", repProgress);

      /* Append bias compensation stats */      
      biasCompStats(preCompPmap, PMAP_TYPE_PRECOMP, tmp);
      strcat(errmsg, tmp);

      rtime = (repLastTime - repStartTime) / 3600.0;
      progress = (float)repProgress / repComplete;
      eta = max(0, rtime * (progress > FTINY ? 1 / progress - 1 : 0));

      sprintf(
         tmp, "%4.2f%% after %.3fu %.3fs %.3fr hours on %s, ETA=%.3fh\n",
         100 * progress, u / 3600, s / 3600, rtime, myhostname(), eta
      );
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
      
      #ifdef SIGCONT
         signal(SIGCONT, pmapPreCompReport);
      #endif
   }



   void pmapDistribReport()
   /* Report photon distribution progress */
   {
      extern char *myhostname();
      float u, s, rtime, eta, progress;
      unsigned t;
      char tmp [512];
      
      #ifdef BSD
         struct rusage rubuf;
      #else
         struct tms tbuf;
         float period;
      #endif

      repLastTime = time(NULL);
      
      #ifdef BSD
         getrusage(RUSAGE_SELF, &rubuf);
         u = rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s = rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
         getrusage(RUSAGE_CHILDREN, &rubuf);
         u += rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s += rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
      #else
         times(&tbuf);
      
         #ifdef _SC_CLK_TCK
            period = 1.0 / sysconf(_SC_CLK_TCK);
         #else
            period = 1.0 / 60;
         #endif
      
         u = (tbuf.tms_utime + tbuf.tms_cutime) * period;
         s = (tbuf.tms_stime + tbuf.tms_cstime) * period;
      #endif

      sprintf(errmsg, "%lu emitted, ", repEmitted);
      
      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (photonMaps [t]) {
            sprintf(tmp, "%lu %s, ", photonMaps [t] -> numPhotons, 
                    pmapName [t]);
            strcat(errmsg, tmp);
         }
      
      rtime = (repLastTime - repStartTime) / 3600.0;
      progress = (float)repProgress / repComplete;
      eta = max(0, rtime * (progress > FTINY ? 1 / progress - 1 : 0));
      
      sprintf(
         tmp, "%4.2f%% after %.3fu %.3fs %.3fr hours on %s, ETA=%.3fh\n", 
         100 * progress, u / 3600, s / 3600, rtime, myhostname(), eta
      );

      strcat(errmsg, tmp);
      eputs(errmsg);
      fflush(stderr);
      
      #ifndef SIGCONT
         signal(SIGCONT, pmapDistribReport);
      #endif
   }
   
#else /* POSIX */

   void pmapPreCompReport()
   /* Report global photon precomputation progress */
   {
      char  tmp [512];
      float rtime, progress, eta;
      
      repLastTime = time(NULL);
      sprintf(errmsg, "%lu precomputed, ", repProgress);

      /* Append bias compensation stats */
      biasCompStats(preCompPmap, PMAP_TYPE_PRECOMP, tmp);
      strcat(errmsg, tmp);      

      rtime = (repLastTime - repStartTime) / 3600.0;
      progress = (float)repProgress / repComplete;
      eta = max(0, rtime * (progress > FTINY ? 1 / progress - 1 : 0));

      sprintf(
         tmp, "%4.2f%% after %5.4f hours, ETA=%.3fh\n", 
         100 * progress, rtime, eta
      );
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
   }



   void pmapDistribReport()
   /* Report photon distribution progress */
   {
      char     tmp [512];
      unsigned t;
      float    rtime, progress, eta;
      
      repLastTime = time(NULL);
      sprintf(errmsg, "%lu emitted, ", repEmitted);

      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (photonMaps [t]) {
            sprintf(tmp, "%lu %s, ", photonMaps [t] -> numPhotons, 
                    pmapName [t]);
            strcat(errmsg, tmp);
         }      

      rtime = (repLastTime - repStartTime) / 3600.0;
      progress = (float)repProgress / repComplete;
      eta = max(0, rtime * (progress > FTINY ? 1 / progress - 1 : 0));

      sprintf(
         tmp, "%4.2f%% after %5.4f hours, ETA=%.3fh\n", 
         100 * progress, rtime, eta
      );
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
   }
#endif



