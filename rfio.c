#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rftime.h"
#include "rfio.h"

struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin)
{
  int i,j,k,l,flag=0,status,msub,ibin,nadd;
  char filename[128],header[256],nfd[32];
  FILE *file;
  struct spectrogram s;
  float *z;
  int nch,j0,j1;
  double freq,samp_rate;
  float length;
  int nchan;
  float s1,s2;

  // Open first file to get number of channels
  sprintf(filename,"%s_%06d.bin",prefix,isub);
	
  // Open file
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("%s does not exist\n",filename);
    s.nsub=0;
    return s;
  }

  // Read header
  status=fread(header,sizeof(char),256,file);
  status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\n",s.nfd0,&s.freq,&s.samp_rate,&length,&nch);

  // Close file
  fclose(file);

  // Compute plotting channel
  if (f0>0.0 && df0>0.0) {
    s.nchan=(int) (df0/s.samp_rate*(float) nch);
    
    j0=(int) ((f0-0.5*df0-s.freq+0.5*s.samp_rate)*(float) nch/s.samp_rate);
    j1=(int) ((f0+0.5*df0-s.freq+0.5*s.samp_rate)*(float) nch/s.samp_rate);
    
    if (j0<0 || j1>nch) {
      fprintf(stderr,"Requested frequency range out of limits\n");
      s.nsub=0;
      s.nchan=0;
      return s;
    }
  } else {
    s.nchan=nch;
    j0=0;
    j1=s.nchan;
  }

  // Number of subints
  s.nsub=nsub/nbin;

  // Allocate
  s.z=(float *) malloc(sizeof(float)*s.nchan*s.nsub);
  z=(float *) malloc(sizeof(float)*nch);
  s.mjd=(double *) malloc(sizeof(double)*s.nsub);
  s.length=(float *) malloc(sizeof(float)*s.nsub);

  // Initialize
  for (j=0;j<s.nchan*s.nsub;j++)
    s.z[j]=0.0;
  for (j=0;j<s.nsub;j++)
    s.mjd[j]=0.0;

  // Loop over files
  for (k=0,i=0,l=0,ibin=0,nadd=0;l<nsub;k++) {
    // Generate filename
    sprintf(filename,"%s_%06d.bin",prefix,k+isub);

    // Open file
    file=fopen(filename,"r");
    if (file==NULL) {
      printf("%s does not exist\n",filename);
	  break;
    }
    printf("opened %s\n",filename);

    // Loop over contents of file
    for (;l<nsub;l++,ibin++) {
      // Read header
      status=fread(header,sizeof(char),256,file);

      if (status==0)
	break;
      status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\n",nfd,&freq,&samp_rate,&length,&nchan);
      s.mjd[i]+=nfd2mjd(nfd)+0.5*length/86400.0;
      s.length[i]+=length;
      nadd++;

      // Read buffer
      status=fread(z,sizeof(float),nch,file);
      if (status==0)
	break;
      
      // Copy
      for (j=0;j<s.nchan;j++) 
	s.z[i+s.nsub*j]+=z[j+j0];

      // Increment
      if (l%nbin==nbin-1) {
	// Scale
	s.mjd[i]/=(float) nadd;

	for (j=0;j<s.nchan;j++) 
	  s.z[i+s.nsub*j]/=(float) nadd;

	ibin=0;
	nadd=0;
	i++;
      }
    }

    // Close file
    fclose(file);
  }

  // Scale last subint
  s.mjd[i]/=(float) nadd;

  for (j=0;j<s.nchan;j++) 
    s.z[i+s.nsub*j]/=(float) nadd;

  // Swap frequency range
  if (f0>0.0 && df0>0.0) {
    s.freq=f0;
    s.samp_rate=df0;
  }

  // Free 
  free(z);

  return s;
}

void write_spectrogram(struct spectrogram s,char *prefix)
{
  int i,j;
  FILE *file;
  char header[256]="",filename[256],nfd[32];
  float *z;
  double mjd;

  // Allocate
  z=(float *) malloc(sizeof(float)*s.nchan);

  // Generate filename
  sprintf(filename,"%s_%06d.bin",prefix,0);

  // Open file
  file=fopen(filename,"w");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
    // Date
    mjd=s.mjd[i]-0.5*s.length[i]/86400.0;
    mjd2nfd(mjd,nfd);

    // Generate header
    sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nEND\n",nfd,s.freq,s.samp_rate,s.length[i],s.nchan);

    // Copy buffer
    for (j=0;j<s.nchan;j++) 
      z[j]=s.z[i+s.nsub*j];

    // Dump contents
    fwrite(header,sizeof(char),256,file);
    fwrite(z,sizeof(float),s.nchan,file);
  }

  // Close file
  fclose(file);

  // Free
  free(z);

  return;
}
