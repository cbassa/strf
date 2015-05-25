#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>
#include "rfio.h"

#define LIM 128

void dec2sex(double x,char *s,int f,int len);

void usage(void)
{
  printf("rfplot: plot RF observations\n\n");
  printf("-p <path>    Path to file prefix /a/b/c_??????.bin\n");
  printf("-o <file>    Output file name [test_000000.bin]\n");
  printf("-s <start>   Number of starting subintegration [0]\n");
  printf("-l <length>  Number of subintegrations to plot [3600]\n");
  printf("-b <nbin>    Number of subintegrations to bin [1]\n");
  printf("-f <freq>    Frequency to zoom into (Hz)\n");
  printf("-w <bw>      Bandwidth to zoom into (Hz)\n");
  printf("-h           This help\n");

  return;
}

int main(int argc,char *argv[])
{
  struct spectrogram s;
  char path[128],outfile[128]="test";
  int arg=0,nsub=3600,nbin=1,isub=0;
  double f0=0.0,df0=0.0;

  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:o:f:w:s:l:b:h"))!=-1) {
      switch (arg) {
	
      case 'p':
	strcpy(path,optarg);
	break;
	
      case 'o':
	strcpy(outfile,optarg);
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
	break;

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

  // Write data
  write_spectrogram(s,outfile);

  // Free
  //  free(s.z);
  //  free(s.mjd);
  //  free(s.length);

  return 0;
}


// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={"::",",,","hms","  "};

  sign=(x<0 ? '-' : ' ');
  x=3600.*fabs(x);

  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;

  if (len==7) sprintf(s,"%c%02i%c%02i%c%07.4f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==6) sprintf(s,"%c%02i%c%02i%c%06.3f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==5) sprintf(s,"%c%02i%c%02i%c%05.2f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==4) sprintf(s,"%c%02i%c%02i%c%04.1f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==2) sprintf(s,"%c%02i%c%02i%c%02i%c",sign,(int) deg,form[f][0],(int) min,form[f][1],(int) floor(sec),form[f][2]);

  return;
}

