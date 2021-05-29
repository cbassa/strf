#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "rftime.h"
#include "rfio.h"

#define LIM 128
#define NMAX 64



void filter(struct spectrogram s,int site_id,float sigma,char *filename,int graves)
{
  int i,j,k,l;
  float s1,s2,avg,std,dz;
  FILE *file;
  double f;
  int *mask;
  float *sig;

  mask=(int *) malloc(sizeof(int)*s.nchan);
  sig=(float *) malloc(sizeof(float)*s.nchan);
  
  // Open file
  file=fopen(filename,"a");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
    // Set mask
    for (j=0;j<s.nchan;j++)
      mask[j]=1;

    // Iterate to remove outliers
    for (k=0;k<10;k++) {

      // Find average
      for (j=0,s1=s2=0.0;j<s.nchan;j++) {
	if (mask[j]==1) {
	  s1+=s.z[i+s.nsub*j];
	  s2+=1.0;
	}
      }
      avg=s1/s2;
      
      // Find standard deviation
      for (j=0,s1=s2=0.0;j<s.nchan;j++) {
	if (mask[j]==1) {
	  dz=s.z[i+s.nsub*j]-avg;
	  s1+=dz*dz;
	  s2+=1.0;
	}
      }
      std=sqrt(s1/s2);

      // Update mask
      for (j=0,l=0;j<s.nchan;j++) {
	if (fabs(s.z[i+s.nsub*j]-avg)>sigma*std) {
	  mask[j]=0;
	  l++;
	}
      }
    }
       // Reset mask
    for (j=0;j<s.nchan;j++) {
      sig[j]=(s.z[i+s.nsub*j]-avg)/std;
      if (sig[j]>sigma) 
	mask[j]=1;
      else
	mask[j]=0;
    }    

    // Find maximum when points are adjacent
    for (j=0;j<s.nchan-1;j++) {
      if (mask[j]==1 && mask[j+1]==1) {
	if (s.z[i+s.nsub*j]<s.z[i+s.nsub*(j+1)])
	  mask[j]=0;
      }
    }
    for (j=s.nchan-2;j>=0;j--) {
      if (mask[j]==1 && mask[j-1]==1) {
	if (s.z[i+s.nsub*j]<s.z[i+s.nsub*(j-1)])
	  mask[j]=0;
      }
    }

    // Mark points
    for (j=0;j<s.nchan;j++) {
      if (mask[j]==1) {
	f=s.freq-0.5*s.samp_rate+(double) j*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0) {
	  if (graves==0)
	    fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sig[j],site_id);
	  else
	    fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,sig[j],site_id);
	}
      }
    }
  } 

  fclose(file);

  free(mask);
  free(sig);

  return;
}

void usage(void)
{
  printf("rffind: Find signals RF observations\n\n");
  printf("-p <path>    Input path to file /a/b/c_??????.bin\n");
  printf("-s <start>   Number of starting subintegration [0]\n");
  printf("-l <length>  Number of subintegrations to plot [3600]\n");
  printf("-f <freq>    Frequency to zoom into (Hz)\n");
  printf("-w <bw>      Bandwidth to zoom into (Hz)\n");
  printf("-o <offset>  Frequency offset to apply\n");
  printf("-C <site>    Site ID\n");
  printf("-g           GRAVES data\n");
  printf("-S           Sigma limit [default: 5.0]\n");
  printf("-h           This help\n");
}

int main(int argc,char *argv[])
{
  int i,j,k,l,j0,j1,m=2,n;
  struct spectrogram s;
  char path[128];
  int isub=0,nsub=0;
  char *env;
  int site_id=0,graves=0;
  float avg,std;
  int arg=0;
  float sigma=5.0;
  FILE *file;
  double f,f0=0.0,df0=0.0;
  char filename[128]="find.dat";

  // Get site
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    site_id=atoi(env);
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }

  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:f:w:s:l:hC:o:S:g"))!=-1) {
      switch (arg) {
	
      case 'p':
	strcpy(path,optarg);
	break;
	
      case 's':
	isub=atoi(optarg);
	break;

      case 'C':
	site_id=atoi(optarg);
	break;
	
      case 'l':
	nsub=atoi(optarg);
	break;
	
      case 'o':
	strcpy(filename,optarg);
	break;

      case 'f':
	f0=(double) atof(optarg);
	break;

      case 'g':
	graves=1;
	break;
	
      case 'S':
	sigma=atof(optarg);
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

  if (nsub==0) {
    // Read data
    for (i=isub;;i++) {
      s=read_spectrogram(path,i,nsub,f0,df0,1,0.0);

      // Exit on emtpy file
      if (s.nsub==0)
	break;
      
      printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);
    
      // Filter
      filter(s,site_id,sigma,filename,graves);

      // Free
      free_spectrogram(s);
    }
  } else {
    // Read data
    s=read_spectrogram(path,isub,nsub,f0,df0,1,0.0);

    // Exit on emtpy file
    if (s.nsub>0) {
      printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);
      
      // Filter
      filter(s,site_id,sigma,filename,graves);
    }

    // Free
    free_spectrogram(s);
  }

  return 0;
}

