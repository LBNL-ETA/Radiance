#ifndef lint
static const char RCSid[] = "$Id: pmapdump.c,v 2.18 2021/02/18 17:08:50 rschregle Exp $";
#endif

/* 
   ======================================================================
   Dump photon maps as RADIANCE scene description or ASCII point list
   to stdout

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme FARESYS") 
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   (c) Tokyo University of Science,
       supported by the JSPS Grants-in-Aid for Scientific Research 
       (KAKENHI JP19KK0115, "Three-Dimensional Light Flow")   
   ======================================================================
   
   $Id: pmapdump.c,v 2.18 2021/02/18 17:08:50 rschregle Exp $
*/



#include "pmap.h"
#include "pmapio.h"
#include "rtio.h"
#include "resolu.h"
#include "random.h"
#include "math.h"


/* Defaults */
/*    Sphere radius as fraction of avg. intersphere dist */
/*    Relative scale for sphere radius (fudge factor) */
/*    Number of spheres */
#define RADCOEFF 0.05
#define RADSCALE 1.0
#define NSPHERES 10000

/* Format for optional ASCII output as XYZ RGB points */
#define POINTFMT "%g\t%g\t%g\t%g\t%g\t%g\n"

/* RADIANCE material and object defs for each photon type */
typedef struct {
   char *mat, *obj;
} RadianceDef;

const RadianceDef radDefs [] = {
   {  "void glow mat.global\n0\n0\n4 %g %g %g 0\n",
      "mat.global sphere obj.global\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void glow mat.pglobal\n0\n0\n4 %g %g %g 0\n",
      "mat.pglobal sphere obj.pglobal\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void glow mat.caustic\n0\n0\n4 %g %g %g 0\n",
      "mat.caustic sphere obj.caustic\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void glow mat.volume\n0\n0\n4 %g %g %g 0\n",
      "mat.volume sphere obj.volume\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void glow mat.direct\n0\n0\n4 %g %g %g 0\n",
      "mat.direct sphere obj.direct\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void glow mat.contrib\n0\n0\n4 %g %g %g 0\n",
      "mat.contrib sphere obj.contrib\n0\n0\n4 %g %g %g %g\n"
   }
};


/* Default colour codes are as follows:   global         = blue
                                          precomp global = cyan
                                          caustic        = red
                                          volume         = green
                                          direct         = magenta 
                                          contrib        = yellow */
const COLOR colDefs [] = {
   {0.25, 0.25, 2}, {0.1, 1, 1}, {1, 0.1, 0.1}, 
   {0.1, 1, 0.1}, {1, 0.1, 1}, {1, 1, 0.1}
};


static int setBool(char *str, unsigned pos, unsigned *var)
{
   switch ((str) [pos]) {
      case '\0': 
         *var = !*var; 
         break;
      case 'y': case 'Y': case 't': case 'T': case '+': case '1': 
         *var = 1; 
         break;
      case 'n': case 'N': case 'f': case 'F': case '-': case '0': 
         *var = 0; 
         break;
      default: 
         return 0;
   }
   
   return 1;
}


int main (int argc, char** argv)
{
   char           format [MAXFMTLEN];
   RREAL          rad, radScale = RADSCALE, extent, dumpRatio;
   unsigned       arg, j, ptype, dim, fluxCol = 0, points = 0;
   long           numSpheres = NSPHERES;
   COLOR          col = {0, 0, 0};
   FILE           *pmapFile;
   PhotonMap      pm;
   PhotonPrimary  pri;
   Photon         p;
#ifdef PMAP_OOC
   char           leafFname [1024];
#endif

   if (argc < 2) {
      puts("Dump photon maps as RADIANCE scene description "
           "or ASCII point list\n");
      printf("Usage: %s "
             "[-a] [-r radscale1] [-n num1] "
             "[-f | -c rcol1 gcol1 bcol1] pmap1 "
             "[-a] [-r radscale2] [-n num2] "
             "[-f | -c rcol2 gcol2 bcol2] pmap2 "
             "...\n", argv [0]);
      return 1;
   }
   
   for (arg = 1; arg < argc; arg++) {
      /* Parse options */
      if (argv [arg][0] == '-') {
         switch (argv [arg][1]) {
            case 'a':
               if (!setBool(argv [arg], 2, &points))
                  error(USER, "invalid option syntax at -a");
               break;
            case 'r':
               if ((radScale = atof(argv [++arg])) <= 0)
                  error(USER, "invalid radius scale");
               break;
               
            case 'n':
               if ((numSpheres = parseMultiplier(argv [++arg])) <= 0)
                  error(USER, "invalid number of points/spheres");
               break;
               
            case 'c':
               if (fluxCol)
                  error(USER, "-f and -c are mutually exclusive");
               
               if (badarg(argc - arg - 1, &argv [arg + 1], "fff"))
                  error(USER, "invalid RGB colour");
                                 
               for (j = 0; j < 3; j++)
                  col [j] = atof(argv [++arg]);
               break;
               
            case 'f':
               if (intens(col) > 0)
                  error(USER, "-f and -c are mutually exclusive");
                  
               if (!setBool(argv [arg], 2, &fluxCol))
                  error(USER, "invalid option syntax at -f");
               break;
               
            default:
               sprintf(errmsg, "unknown option %s", argv [arg]);
               error(USER, errmsg);
               return -1;
         }
         
         continue;
      }

      /* Open next photon map file */
      if (!(pmapFile = fopen(argv [arg], "rb"))) {
         sprintf(errmsg, "can't open %s", argv [arg]);
         error(SYSTEM, errmsg);
      }

      /* Get format string */
      strcpy(format, PMAP_FORMAT_GLOB);
      if (checkheader(pmapFile, format, NULL) != 1) {
         sprintf(errmsg, "photon map file %s has unknown format %s", 
                 argv [arg], format);
         error(USER, errmsg);
      }

      /* Identify photon map type from format string */
      for (ptype = 0; 
           ptype < NUM_PMAP_TYPES && strcmp(pmapFormat [ptype], format);
           ptype++);

      if (!validPmapType(ptype)) {
         sprintf(errmsg, "file %s contains an unknown photon map type", 
                argv [arg]);
         error(USER, errmsg);
      }

      /* Get file format version and check for compatibility */
      if (strcmp(getstr(format, pmapFile), PMAP_FILEVER))      
         error(USER, "incompatible photon map file format");

      if (!points) {
         /* Dump command line as comment */
         fputs("# ", stdout);
         printargs(argc, argv, stdout);
         fputc('\n', stdout);
      }
         
      /* Set point/sphere colour if independent of photon flux, 
         output RADIANCE material def if required */
      if (!fluxCol) {
         if (intens(col) <= 0)
            copycolor(col, colDefs [ptype]);
         if (!points) {
            printf(radDefs [ptype].mat, col [0], col [1], col [2]);
            fputc('\n', stdout);
         }
      }

      /* Get number of photons as fixed size, which possibly results in
       * padding of MSB with 0 on some platforms.  Unlike sizeof() however,
       * this ensures portability since this value may span 32 or 64 bits
       * depending on platform.  */
      pm.numPhotons = getint(PMAP_LONGSIZE, pmapFile);      

      /* Skip avg photon flux */ 
      for (j = 0; j < 3; j++) 
         getflt(pmapFile);

      /* Get distribution extent (min & max photon positions) */
      for (j = 0; j < 3; j++) {
         pm.minPos [j] = getflt(pmapFile);
         pm.maxPos [j] = getflt(pmapFile);
      }
      
      /* Skip centre of gravity, and avg photon dist to it */
      for (j = 0; j < 4; j++)
         getflt(pmapFile);
      
      /* Sphere radius based on avg intersphere dist depending on
         whether the distribution occupies a 1D line (!), a 2D plane, 
         or 3D volume (= sphere distrib density ^-1/d, where d is the
         dimensionality of the distribution) */
      for (j = 0, extent = 1.0, dim = 0; j < 3; j++) {
         rad = pm.maxPos [j] - pm.minPos [j];
         
         if (rad > FTINY) {
            dim++;
            extent *= rad;
         }
      }

      rad = radScale * RADCOEFF * pow(extent / numSpheres, 1./dim);
      
      /* Photon dump probability to satisfy target sphere count */
      dumpRatio = min(1, (float)numSpheres / pm.numPhotons);
      
      /* Skip primary rays */
      pm.numPrimary = getint(sizeof(pm.numPrimary), pmapFile);
      while (pm.numPrimary-- > 0) {
         /* Skip source index & incident dir */
         getint(sizeof(pri.srcIdx), pmapFile);
#ifdef PMAP_PRIMARYDIR
         /* Skip primary incident dir */
         getint(sizeof(pri.dir), pmapFile);         
#endif
#ifdef PMAP_PRIMARYPOS
         /* Skip primary hitpoint */
         for (j = 0; j < 3; j++)
            getflt(pmapFile);
#endif
      }

#ifdef PMAP_OOC
      /* Open leaf file with filename derived from pmap, replace pmapFile
       * (which is currently the node file) */
      strncpy(leafFname, argv [arg], sizeof(leafFname) - 1);
      strncat(leafFname, PMAP_OOC_LEAFSUFFIX, sizeof(leafFname) - 1);
      fclose(pmapFile);
      if (!(pmapFile = fopen(leafFname, "rb"))) {
         sprintf(errmsg, "cannot open leaf file %s", leafFname);
         error(SYSTEM, errmsg);
      }
#endif
            
      /* Read photons */
      while (pm.numPhotons-- > 0) {
#ifdef PMAP_OOC
         /* Get entire photon record from ooC octree leaf file
            !!! OOC PMAP FILES CURRENTLY DON'T USE PORTABLE I/O !!! */
         if (!fread(&p, sizeof(p), 1, pmapFile)) {
            sprintf(errmsg, "error reading OOC leaf file %s", leafFname);
            error(SYSTEM, errmsg);
         }
#else /* kd-tree */
         /* Get photon position */
         for (j = 0; j < 3; j++) 
            p.pos [j] = getflt(pmapFile);

         /* Get photon normal (currently not used) */
         for (j = 0; j < 3; j++) 
            p.norm [j] = getint(1, pmapFile);

         /* Get photon flux */
   #ifdef PMAP_FLOAT_FLUX
         for (j = 0; j < 3; j++) 
            p.flux [j] = getflt(pmapFile);
   #else
         for (j = 0; j < 4; j++) 
            p.flux [j] = getint(1, pmapFile);
   #endif
   
         

         /* Skip primary ray index */
         getint(sizeof(p.primary), pmapFile);

         /* Skip flags */
         getint(sizeof(p.flags), pmapFile);
#endif

         /* Dump photon probabilistically acc. to target sphere count */
         if (frandom() <= dumpRatio) {
            if (fluxCol) {
               /* Get photon flux */
               getPhotonFlux(&p, col);
               /* Scale by dumpRatio for energy conservation */
               scalecolor(col, 1.0 / dumpRatio);
            }
            
            if (!points) {
               if (fluxCol) {
                  /* Dump material def if variable (depends on flux) */
                  printf(radDefs [ptype].mat, col [0], col [1], col [2]);
                  fputc('\n', stdout);
               }
               printf(radDefs [ptype].obj, p.pos [0], p.pos [1], p.pos [2],
                      rad);
               fputc('\n', stdout);
            }
            else /* Dump as XYZ RGB point */
               printf(POINTFMT, p.pos [0], p.pos [1], p.pos [2],
                      col [0], col [1] ,col [2]);
         }
         
         if (ferror(pmapFile) || feof(pmapFile)) {
            sprintf(errmsg, "error reading %s", argv [arg]);
            error(USER, errmsg);
         }
      }
      
      fclose(pmapFile);
      
      /* Reset defaults for next dump */
      radScale = RADSCALE;
      numSpheres = NSPHERES;
      col [0] = col [1] = col [2] = 0;
      fluxCol = points = 0;
   }
   
   return 0;
}
