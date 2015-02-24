/* 
   ==================================================================
   Photon map generator
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: mkpmap.c,v 2.1 2015/02/24 19:39:26 greg Exp $    
*/



#include "pmap.h"
#include "pmapmat.h"
#include "pmapcontrib.h"
#include "pmaprand.h"
#include "paths.h"
#include "ambient.h"
#include "resolu.h"
#include "source.h"
#include <string.h>
#include <sys/stat.h>



extern char VersionID [];



char* progname;                         /* argv[0] */
int  dimlist [MAXDIM];                  /* sampling dimensions */
int  ndims = 0;                         /* number of sampling dimenshunns */
char* octname = NULL;                   /* octree name */
CUBE thescene;                          /* scene top-level octree */
OBJECT nsceneobjs;                      /* number of objects in scene */
double srcsizerat = 0.01;               /* source partition size ratio */
int backvis = 1;                        /* back face visibility */
int clobber = 0;                        /* overwrite output */
COLOR cextinction = BLKCOLOR;           /* global extinction coefficient */
COLOR salbedo = BLKCOLOR;               /* global scattering albedo */
double seccg = 0;                       /* global scattering eccentricity */
int ambincl = -1;                       /* photon port flag */
char *amblist [AMBLLEN + 1];            /* photon port list */
char *diagFile = NULL;                  /* diagnostics output file */



/* Dummies for linkage */

COLOR ambval = BLKCOLOR;
double shadthresh = .05, ambacc = 0.2, shadcert = .5, minweight = 5e-3, 
       ssampdist = 0, dstrsrc = 0.0, specthresh = 0.15, specjitter = 1.0,
       avgrefl = 0.5;    
int ambvwt = 0, ambssamp = 0, ambres = 32, ambounce = 0, directrelay = 1,
    directvis = 1, samplendx, do_irrad = 0, ambdiv = 128, vspretest = 512,
    maxdepth = 6, contrib = 0;
char *shm_boundary = NULL, *ambfile = NULL, *RCCONTEXT = NULL;
void (*trace)() = NULL, (*addobjnotify [])() = {ambnotify, NULL};


  
void printdefaults()
/* print default values to stdout */
{
   puts("-apg file nPhotons\t\t# global photon map");
   puts("-apc file nPhotons\t\t# caustic photon map");          
   puts("-apd file nPhotons\t\t# direct photon map");
   puts("-app file nPhotons bwidth\t# precomputed global photon map");
#if 0
   /* Hide this option as most likely useless and confusing to user */
   puts("-appb file nPhotons minBw maxBw\t# precomp bias comp global pmap");
#endif
   puts("-apv file nPhotons\t\t# volume photon map");
   puts("-apC file nPhotons\t\t# contribution photon map");
   
   printf("-apD %f\t\t\t# predistribution factor\n", preDistrib);
   printf("-apM %d\t\t\t\t# max predistrib passes\n", maxPreDistrib);
   printf("-apm %ld\t\t\t# max photon bounces\n", photonMaxBounce);                            
   puts("-apo mod\t\t\t# photon port modifier");
   puts("-apO file\t\t\t# photon port file");
   printf("-apP %f\t\t\t# precomputation factor\n", finalGather);
   printf("-apr %d\t\t\t\t# random seed\n", randSeed);
   puts("-aps mod\t\t\t# antimatter sensor modifier");
   puts("-apS file\t\t\t# antimatter sensor file");

   printf(backvis ? "-bv+\t\t\t\t# back face visibility on\n"
                  : "-bv-\t\t\t\t# back face visibility off\n");
   printf("-dp  %.1f\t\t\t# PDF samples / sr\n", pdfSamples);
   printf("-ds  %f\t\t\t# source partition size ratio\n", srcsizerat);
   printf("-e   %s\t\t\t# diagnostics output file\n", diagFile);
   printf(clobber ? "-fo+\t\t\t\t# force overwrite"
                  : "-fo-\t\t\t\t# do not overwrite\n");
   printf("-i   %-9ld\t\t\t# photon heap size increment\n", 
          photonHeapSizeInc);
   printf("-ma  %.2f %.2f %.2f\t\t# scattering albedo\n", 
          colval(salbedo,RED), colval(salbedo,GRN), colval(salbedo,BLU));
   printf("-me  %.2e %.2e %.2e\t# extinction coefficient\n", 
          colval(cextinction,RED), colval(cextinction,GRN), 
          colval(cextinction,BLU));          
   printf("-mg  %.2f\t\t\t# scattering eccentricity\n", seccg);
   printf("-t   %-9d\t\t\t# time between reports\n", photonRepTime);

#ifdef PMAP_ROI
   /* Ziss option for ze egg-spurtz only! */
   printf("-api \t%.0e %.0e %.0e\n\t%.0e %.0e %.0e\t# region of interest\n", 
          pmapROI [0], pmapROI [1], pmapROI [2], pmapROI [3], 
          pmapROI [4], pmapROI [5]);
#endif   
}




int main (int argc, char* argv [])
{
   #define check(ol, al) if (argv [i][ol] || \
                             badarg(argc - i - 1,argv + i + 1, al)) \
                            goto badopt
                            
   #define bool(olen, var) switch (argv [i][olen]) { \
                             case '\0': var = !var; break; \
                             case 'y': case 'Y': case 't': case 'T': \
                             case '+': case '1': var = 1; break; \
                             case 'n': case 'N': case 'f': case 'F': \
                             case '-': case '0': var = 0; break; \
                             default: goto badopt; \
                          }   

   int loadflags = IO_CHECK | IO_SCENE | IO_TREE | IO_BOUNDS, rval, i;
   char **portLp = NULL, **sensLp = photonSensorList;
   struct stat pmstat;

   /* Global program name */
   progname = fixargv0(argv [0]);
   /* Initialize object types */
   initotypes();
   
   /* Parse options */
   for (i = 1; i < argc; i++) {
      /* Eggs-pand arguments */
      while ((rval = expandarg(&argc, &argv, i)))
         if (rval < 0) {
            sprintf(errmsg, "cannot eggs-pand '%s'", argv [i]);
            error(SYSTEM, errmsg);
         }
         
      if (argv[i] == NULL) 
         break;
         
      if (!strcmp(argv [i], "-version")) {
         puts(VersionID);
         quit(0);
      }
      
      if (!strcmp(argv [i], "-defaults") || !strcmp(argv [i], "-help")) {
         printdefaults();
         quit(0);
      }
      
      /* Get octree */
      if (i == argc - 1) {
         octname = argv [i];
         break;
      }
            
      switch (argv [i][1]) {
         case 'a': 
            if (!strcmp(argv [i] + 2, "pg")) {
               /* Global photon map */
               check(4, "ss");
               globalPmapParams.fileName = argv [++i];
               globalPmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!globalPmapParams.distribTarget) 
                  goto badopt;                         
               globalPmapParams.minGather = globalPmapParams.maxGather = 0;
            }
                              
            else if (!strcmp(argv [i] + 2, "pm")) {
               /* Max photon bounces */
               check(4, "i");
               photonMaxBounce = atol(argv [++i]);
               if (!photonMaxBounce) 
                  goto badopt;
            }
            
            else if (!strcmp(argv [i] + 2, "pp")) {
               /* Precomputed global photon map */
               check(4, "ssi");
               preCompPmapParams.fileName = argv [++i];
               preCompPmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!preCompPmapParams.distribTarget) 
                  goto badopt;
               preCompPmapParams.minGather = preCompPmapParams.maxGather = 
                  atoi(argv [++i]);
               if (!preCompPmapParams.maxGather) 
                  goto badopt;
            }
            
#if 0
            else if (!strcmp(argv [i] + 2, "ppb")) {
               /* Precomputed global photon map + bias comp. */
               check(5, "ssii");
               preCompPmapParams.fileName = argv [++i];
               preCompPmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!preCompPmapParams.distribTarget) 
                  goto badopt;                      
               preCompPmapParams.minGather = atoi(argv [++i]);
               preCompPmapParams.maxGather = atoi(argv [++i]);
               if (!preCompPmapParams.minGather ||
                   preCompPmapParams.minGather >= 
                   preCompPmapParams.maxGather) 
                  goto badopt;
            }
#endif
            
            else if (!strcmp(argv [i] + 2, "pc")) {
               /* Caustic photon map */
               check(4, "ss");
               causticPmapParams.fileName = argv [++i];
               causticPmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!causticPmapParams.distribTarget) 
                  goto badopt;
            }
            
            else if (!strcmp(argv [i] + 2, "pv")) {
               /* Volume photon map */
               check(4, "ss");
               volumePmapParams.fileName = argv [++i];
               volumePmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!volumePmapParams.distribTarget) 
                  goto badopt;                      
            }
            
            else if (!strcmp(argv [i] + 2, "pd")) {
               /* Direct photon map */
               check(4, "ss");
               directPmapParams.fileName = argv [++i];
               directPmapParams.distribTarget = 
                  parseMultiplier(argv [++i]);
               if (!directPmapParams.distribTarget) 
                  goto badopt;
            }
            
            else if (!strcmp(argv [i] + 2, "pC")) {
               /* Light source contribution photon map */
               check(4, "ss");
               contribPmapParams.fileName = argv [++i];
               contribPmapParams.distribTarget =
                  parseMultiplier(argv [++i]);
               if (!contribPmapParams.distribTarget)
                  goto badopt;
            }
            
            else if (!strcmp(argv [i] + 2, "pD")) {
               /* Predistribution factor */
               check(4, "f");
               preDistrib = atof(argv [++i]);
               if (preDistrib <= 0)
                  error(USER, "predistribution factor must be > 0");
            }   

#ifdef PMAP_ROI
            /* Region of interest; ziss option for ze egg-spurtz only! */
            else if (!strcmp(argv [i] + 2, "pi")) {
               int j;
               check(4, "ffffff");
               for (j = 0; j < 6; j++)
                  pmapROI [j] = atof(argv [++i]);
            }               
#endif
             
            else if (!strcmp(argv [i] + 2, "pP")) {
               /* Global photon precomputation factor */
               check(4, "f");
               finalGather = atof(argv [++i]);
               if (finalGather <= 0 || finalGather > 1)
                  error(USER, "global photon precomputation factor "
                        "must be in range ]0, 1]");
            }                  
            
            else if (!strcmp(argv [i] + 2, "po") || 
                     !strcmp(argv [i] + 2, "pO")) {
               /* Photon port */
               check(4, "s");
               
               if (ambincl != 1) {
                  ambincl = 1;
                  portLp = amblist;
               }
               
               if (argv[i][3] == 'O') {	
                  /* Get port modifiers file */
                  rval = wordfile(portLp, getpath(argv [++i],
                                  getrlibpath(), R_OK));
                  if (rval < 0) {
                      sprintf(errmsg, "cannot open photon port file %s",
                              argv [i]);
                      error(SYSTEM, errmsg);
                  }
                  
                  portLp += rval;
               } 
               
               else {
                  /* Append modifier to port list */
                  *portLp++ = argv [++i];
                  *portLp = NULL;
               }
            }
            
            else if (!strcmp(argv [i] + 2, "pr")) {
               /* Random seed */
               check(4, "i");
               randSeed = atoi(argv [++i]);
            }                   

            else if (!strcmp(argv [i] + 2, "ps") || 
                     !strcmp(argv [i] + 2, "pS")) {
               /* Antimatter sensor */
               check(4, "s");
               
               if (argv[i][3] == 'S') {	
                  /* Get sensor modifiers from file */
                  rval = wordfile(sensLp, getpath(argv [++i],
                                  getrlibpath(), R_OK));
                  if (rval < 0) {
                      sprintf(errmsg, "cannot open antimatter sensor file %s",
                              argv [i]);
                      error(SYSTEM, errmsg);
                  }
                  
                  sensLp += rval;
               } 
               
               else {
                  /* Append modifier to sensor list */
                  *sensLp++ = argv [++i];
                  *sensLp = NULL;
               }
            }

            else goto badopt;                   
            break;
                   
         case 'b': 
            if (argv [i][2] == 'v') {
               /* Back face visibility */
               bool(3, backvis);
            }
                   
            else goto badopt;
            break;
                   
         case 'd': /* Direct */
            switch (argv [i][2]) {
               case 'p': /* PDF samples */
                  check(3, "f");
                  pdfSamples = atof(argv [++i]);
                  break;
                  
               case 's': /* Source partition size ratio */
                  check(3, "f");
                  srcsizerat = atof(argv [++i]);
                  break;
                  
               default: goto badopt;
            }                   
            break;
                   
         case 'e': /* Diagnostics file */
            check(2, "s");
            diagFile = argv [++i];
            break;
                  
         case 'f': 
            if (argv [i][2] == 'o') {
               /* Force overwrite */
               bool(3, clobber);
            }
                   
            else goto badopt;
            break; 
            
         case 'i': /* Photon heap size increment */
            check(2, "i");
            photonHeapSizeInc = atol(argv [++i]);
            break;
                   
         case 'm': /* Medium */
            switch (argv[i][2]) {
               case 'e':	/* Eggs-tinction coefficient */
                  check(3, "fff");
                  setcolor(cextinction, atof(argv [i + 1]),
                           atof(argv [i + 2]), atof(argv [i + 3]));
                  i += 3;
                  break;
                                
               case 'a':	/* Albedo */
                  check(3, "fff");
                  setcolor(salbedo, atof(argv [i + 1]), 
                           atof(argv [i + 2]), atof(argv [i + 3]));
                  i += 3;
                  break;
                                
               case 'g':	/* Scattering eccentricity */
                  check(3, "f");
                  seccg = atof(argv [++i]);
                  break;
                                
               default: goto badopt;
            }                   
            break;
                   
         case 't': /* Timer */
            check(2, "i");
            photonRepTime = atoi(argv [++i]);
            break;
                   
         default: goto badopt;
      }
   }
   
   /* Open diagnostics file */
   if (diagFile) {
      if (!freopen(diagFile, "a", stderr)) quit(2);
      fprintf(stderr, "**************\n*** PID %5d: ", getpid());
      printargs(argc, argv, stderr);
      putc('\n', stderr);
      fflush(stderr);
   }
   
#ifdef NICE
   /* Lower priority */
   nice(NICE);
#endif

   if (octname == NULL) 
      error(USER, "missing octree argument");
      
   /* Allocate photon maps and set parameters */
   for (i = 0; i < NUM_PMAP_TYPES; i++) {
      setPmapParam(photonMaps + i, pmapParams + i);
      
      /* Don't overwrite existing photon map unless clobbering enabled */
      if (photonMaps [i] && !stat(photonMaps [i] -> fileName, &pmstat) &&
          !clobber) {
         sprintf(errmsg, "photon map file %s exists, not overwritten",
                 photonMaps [i] -> fileName);
         error(USER, errmsg);
      }
   }
       
   for (i = 0; i < NUM_PMAP_TYPES && !photonMaps [i]; i++);   
   if (i >= NUM_PMAP_TYPES)
      error(USER, "no photon maps specified");
   
   readoct(octname, loadflags, &thescene, NULL);
   nsceneobjs = nobjects;
   
   /* Get sources */
   marksources();

   /* Do forward pass and build photon maps */
   if (contribPmap)
      /* Just build contrib pmap, ignore others */
      distribPhotonContrib(contribPmap);
   else
      distribPhotons(photonMaps);
   
   /* Save photon maps; no idea why GCC needs an explicit cast here... */
   savePmaps((const PhotonMap**)photonMaps, argc, argv);
   cleanUpPmaps(photonMaps);
   
   quit(0);

badopt:
   sprintf(errmsg, "command line error at '%s'", argv[i]);
   error(USER, errmsg);

   #undef check
   #undef bool
   return 0;
}
