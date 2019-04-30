#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include "rfio.h"
#include "rftime.h"
#include "rftrace.h"

void usage(void)
{
  printf("rfdop m:t:c:l:d:g:s:o:hi:\n\n");
  printf("m    date/time (MJD)\n");
  printf("t    date/time (yyyy-mm-ddThh:mm:ss.sss)\n");
  printf("i    NORAD number\n");
  printf("c    TLE catalog file\n");
  printf("s    site (COSPAR)\n");
  printf("d    time step [default: 10s]\n");
  printf("l    trail length [default: 900s]\n");
  printf("g    enable GRAVES\n");
  printf("H    skip high orbits (<10 revs/day)\n");
  printf("o    output file name\n");
  printf("h    this help\n");
  
  return;
}

int main(int argc, char *argv[])
{
  int i, n, satno=0;
  double mjd0, *mjd, t;
  float length=900.0, dt=10.0;
  int site_id;
  int graves=0,skiphigh=0;
  char *env,tlefile[256],nfd[64],outfname[256]="out.dat";
  int arg=0;

  // Get defaults
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    site_id=atoi(env);
  }

  if (argc>1) {
    while ((arg=getopt(argc,argv,"m:t:c:i:l:d:gs:o:hH"))!=-1) {
      switch (arg) {
      case 't':
	strcpy(nfd,optarg);
	mjd0=nfd2mjd(nfd);
	break;
	
      case 'm':
	mjd0=(double) atof(optarg);
	break;
	
      case 'c':
	strcpy(tlefile,optarg);
	break;

      case 'i':
	satno=atoi(optarg);
	break;
	
      case 'l':
	length=atof(optarg);
	break;
	
      case 'd':
	dt=atof(optarg);
	break;
	
      case 's':
	site_id=atoi(optarg);
	break;
	
      case 'g':
	graves=1;
	break;

      case 'H':
	skiphigh=1;
	break;
	
      case 'o':
	strcpy(outfname,optarg);
	break;

      case 'h':
	usage();
	return 0;
	
      default:
	usage();
	return 0;
      }
    }
  } else {
    usage();
    return 0;
  }
  
  // Generate MJDs
  n = (int) floor(length/dt);
  mjd = (double *) malloc(sizeof(double)*n);
  for (i=0;i<n;i++) 
    mjd[i] = mjd0+(double) i*dt/86400.0;

  // Compute Doppler curves
  compute_doppler(tlefile,mjd,n,site_id,satno,graves,skiphigh,outfname);
  
  // Free
  free(mjd);
    
  return 0;
}
