/* 
   ==================================================================
   Photon map file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: pmapio.c,v 2.1 2015/02/24 19:39:27 greg Exp $    
*/



#include "pmapio.h"
#include "pmapdiag.h"
#include "resolu.h"



void savePhotonMap (const PhotonMap *pmap, const char *fname,
                    PhotonMapType type, int argc, char **argv)
{
   unsigned long i, j;
   const Photon* p;
   FILE* file;

   if (!pmap || !pmap -> heap || !pmap -> heapSize || 
       !validPmapType(type)) {
      error(INTERNAL, "attempt to save empty or invalid photon map");
      return;
   }
      
   if (photonRepTime) {
      sprintf(errmsg, "Saving %s (%ld photons)...\n", 
              fname, pmap -> heapSize);
      eputs(errmsg);
      fflush(stderr);
   }
   
   if (!(file = fopen(fname, "wb"))) {
      sprintf(errmsg, "can't open photon map file %s", fname);
      error(SYSTEM, errmsg);
   }
      
   /* Write header */
   newheader("RADIANCE", file);
   
   /* Write command line */
   printargs(argc, argv, file);
   
   /* Include statistics in info text */
   fprintf(file, "%ld photons @ (%.2e, %.2e, %.2e) avg watts\n"
           "Extent [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f]\n",
           pmap -> heapSize, pmap -> photonFlux [0], 
           pmap -> photonFlux [1], pmap -> photonFlux [2],
           pmap -> minPos [0], pmap -> minPos [1], pmap -> minPos [2],
           pmap -> maxPos [0], pmap -> maxPos [1], pmap -> maxPos [2]);
   if (pmap -> primary)
      fprintf(file, "%d primary rays\n", pmap -> primaryEnd + 1);
   
   /* Write format */
   fputformat((char*)pmapFormat [type], file);
   fprintf(file, "VERSION=%d\n", PMAP_FILEVER);
   
   /* Empty line = end of header */
   putc('\n', file);
   
   /* Write file format version */
   putint(PMAP_FILEVER, sizeof(PMAP_FILEVER), file); 
   
   /* Write number of photons */
   putint(pmap -> heapSize, sizeof(pmap -> heapSize), file);
   
   /* Write average photon flux */
   for (j = 0; j < 3; j++) 
      putflt(pmap -> photonFlux [j], file);

   /* Write max and min photon positions */
   for (j = 0; j < 3; j++) {
      putflt(pmap -> minPos [j], file);
      putflt(pmap -> maxPos [j], file);
   }

   /* Write centre of gravity */
   for (j = 0; j < 3; j++)
      putflt(pmap -> CoG [j], file);
      
   /* Write avg distance to centre of gravity */
   putflt(pmap -> CoGdist, file);
         
   for (i = 0, p = pmap -> heap; i < pmap -> heapSize; i++, p++) {
      /* Write photon attributes */
      for (j = 0; j < 3; j++)
         putflt(p -> pos [j], file);
         
      /* Bytewise dump otherwise we have portability probs */
      for (j = 0; j < 3; j++) 
         putint(p -> norm [j], 1, file);
         
      #ifdef PMAP_FLOAT_FLUX
         for (j = 0; j < 3; j++) 
            putflt(p -> flux [j], file);
      #else
         for (j = 0; j < 4; j++) 
            putint(p -> flux [j], 1, file);
      #endif

      putint(p -> primary, sizeof(p -> primary), file);
      putint(p -> flags, 1, file);
            
      if (ferror(file)) {
         sprintf(errmsg, "error writing photon map file %s", fname);
         error(SYSTEM, errmsg);
      }
   }
   
   /* Write out primary photon rays (or just zero count if none) */
   if (pmap -> primary) {
      /* primaryEnd points to last primary ray in array, so increment for
       * number of entries */
      putint(pmap -> primaryEnd + 1, sizeof(pmap -> primaryEnd), file);
      
      for (i = 0; i <= pmap -> primaryEnd; i++) {
         PhotonPrimary *prim = pmap -> primary + i;      
         
         putint(prim -> srcIdx, sizeof(prim -> srcIdx), file);

         for (j = 0; j < 3; j++)
            putflt(prim -> dir [j], file);
            
         for (j = 0; j < 3; j++)
            putflt(prim -> org [j], file);
         
         if (ferror(file))
            error(SYSTEM, "error writing primary photon rays");
      }
   }
   else putint(0, sizeof(pmap -> primaryEnd), file);

   fclose(file);
}



PhotonMapType loadPhotonMap (PhotonMap *pmap, const char *fname)
{
   Photon* p;
   PhotonMapType ptype = PMAP_TYPE_NONE;
   char format [128];
   unsigned long i, j;
   FILE *file;

   if (!pmap)
      return PMAP_TYPE_NONE;
   
   if ((file = fopen(fname, "rb")) == NULL) {
      sprintf(errmsg, "can't open photon map file %s", fname);
      error(SYSTEM, errmsg);
   }

   /* Get format string */
   strcpy(format, PMAP_FORMAT_GLOB);
   if (checkheader(file, format, NULL) != 1) {
      sprintf(errmsg, "photon map file %s has unknown format %s", 
              fname, format);
      error(USER, errmsg);
   }
   
   /* Identify photon map type from format string */
   for (ptype = 0; 
        strcmp(pmapFormat [ptype], format) && ptype < NUM_PMAP_TYPES;
        ptype++);
      
   if (!validPmapType(ptype)) {
      sprintf(errmsg, "file %s contains an unknown photon map type", fname);
      error(USER, errmsg);
   }

   initPhotonMap(pmap, ptype);
   
   /* Get file format version and check for compatibility */
   if (getint(sizeof(PMAP_FILEVER), file) != PMAP_FILEVER)
      error(USER, "incompatible photon map file format");
      
   /* Get number of photons */
   pmap -> heapSize = pmap -> heapEnd = 
      getint(sizeof(pmap -> heapSize), file);
   pmap -> heap = (Photon *)malloc(pmap -> heapSize * sizeof(Photon));
   if (!pmap -> heap) 
      error(INTERNAL, "can't allocate photon heap from file");
      
   /* Get average photon flux */
   for (j = 0; j < 3; j++) 
      pmap -> photonFlux [j] = getflt(file);
   
   /* Get max and min photon positions */
   for (j = 0; j < 3; j++) {
      pmap -> minPos [j] = getflt(file);
      pmap -> maxPos [j] = getflt(file);
   }
   
   /* Get centre of gravity */
   for (j = 0; j < 3; j++)
      pmap -> CoG [j] = getflt(file);
      
   /* Get avg distance to centre of gravity */
   pmap -> CoGdist = getflt(file);
      
   /* Get photon attributes */   
   for (i = 0, p = pmap -> heap; i < pmap -> heapSize; i++, p++) {
      for (j = 0; j < 3; j++) 
         p -> pos [j] = getflt(file);
         
      /* Bytewise grab otherwise we have portability probs */
      for (j = 0; j < 3; j++) 
         p -> norm [j] = getint(1, file);
         
      #ifdef PMAP_FLOAT_FLUX
         for (j = 0; j < 3; j++) 
            p -> flux [j] = getflt(file);
      #else      
         for (j = 0; j < 4; j++) 
            p -> flux [j] = getint(1, file);
      #endif

      p -> primary = getint(sizeof(p -> primary), file);
      p -> flags = getint(1, file);
      
      if (feof(file)) {
         sprintf(errmsg, "error reading photon map file %s", fname);
         error(SYSTEM, errmsg);
      }
   }
      
   /* Get primary photon rays */
   pmap -> primarySize = getint(sizeof(pmap -> primarySize), file);
   
   if (pmap -> primarySize) {
      pmap -> primaryEnd = pmap -> primarySize - 1;
      pmap -> primary = (PhotonPrimary*)malloc(pmap -> primarySize *
                                               sizeof(PhotonPrimary));
      if (!pmap -> primary)
         error(INTERNAL, "can't allocate primary photon rays");
         
      for (i = 0; i < pmap -> primarySize; i++) {
         PhotonPrimary *prim = pmap -> primary + i;
         
         prim -> srcIdx = getint(sizeof(prim -> srcIdx), file);
         
         for (j = 0; j < 3; j++)
            prim -> dir [j] = getflt(file);
            
         for (j = 0; j < 3; j++)
            prim -> org [j] = getflt(file);
            
         if (feof(file))
            error(SYSTEM, "error reading primary photon rays");
      }
   }
   
   fclose(file);
   
   return ptype;
}
