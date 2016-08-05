#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>
#include "rftime.h"
#include "rfio.h"

void usage(void);

int main(int argc,char *argv[])
{
  int arg=0;
  struct spectrogram s;
  char path[128],filename[32];
  int isub=0,nsub=3600,nbin=1,jsub;
  double f0=0.0,df0=0.0;
  float fmin,fmax,f;
  float zmin=0.0,zmax=0.12,z;
  int ichan,k;

  zmin=-30.0;
  zmax=15.0;
  
  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:hc:C:gm:o:"))!=-1) {
      switch (arg) {
	
      case 'p':
	strcpy(path,optarg);
	break;
	
      case 's':
	isub=atoi(optarg);
	break;
	
      case 'l':
	nsub=atoi(optarg);
	break;
	
      case 'b':
	nbin=atoi(optarg);
	break;
	
      case 'f':
	f0=(double) atof(optarg);
	break;
	
      case 'w':
	df0=(double) atof(optarg);
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

  // Read data
  s=read_spectrogram(path,isub,nsub,f0,df0,nbin,0.0);

  printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

  // Exit on empty data
  if (s.nsub==0)
    return 0;

  // Open plot
  cpgopen("/xs");
  cpgsch(0.8);

  // Frequency axis
  fmin=-0.5*s.samp_rate*1e-6;
  fmax=0.5*s.samp_rate*1e-6;

  cpgenv(fmin,fmax,zmin,zmax,0,0);
  cpglab("Frequency offset (MHz)","Power (dB)"," ");

  for (jsub=0;jsub<nsub;jsub++) {
    for (ichan=0;ichan<s.nchan;ichan++) {
      f=-0.5*s.samp_rate+s.samp_rate*(float) ichan/(float) s.nchan;
      f*=1e-6;
      k=isub+jsub+s.nsub*ichan;
      z=10.0*log10(s.z[k]+1e-20);
      if (ichan==0)
	cpgmove(f,z);
      else
	cpgdraw(f,z);
    }
  }
  cpgend();

  // Free
  free(s.z);
  free(s.zavg);
  free(s.zstd);
  free(s.mjd);

  return 0;
}

void usage(void)
{
  printf("rfspec: plot RF spectra\n\n");
  printf("-p <path>    Input path to file /a/b/c_??????.bin\n");
  printf("-s <start>   Number of starting subintegration [0]\n");
  printf("-l <length>  Number of subintegrations to plot [3600]\n");
  printf("-b <nbin>    Number of subintegrations to bin [1]\n");
  printf("-f <freq>    Frequency to zoom into (Hz)\n");
  printf("-w <bw>      Bandwidth to zoom into (Hz)\n");
  printf("-h           This help\n");

  return;
}
