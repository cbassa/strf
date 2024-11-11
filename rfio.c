#define _XOPEN_SOURCE
#include "rfio.h"
#include "rftime.h"
#include "zscale.h"
#include <dirent.h>
#include <libgen.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


// Get the first and last bin number matching the prefix. There is no
// guarantee that the range is continous and there is no missing bin in
// the range.
int get_bin_range(char *prefix, int *first_bin, int *last_bin) {
  DIR *dp;
  struct dirent *ep;

  // Both, dirname and basename, don't handle static strings and can modifiy
  // their argument. So make a copy and use a pointer to the result.
  char dir[256];
  char base[256];
  strncpy(dir, prefix, 256);
  strncpy(base, prefix, 256);
  char *dir_p = dirname(dir);
  char *base_p = strcat(basename(base), "_");

  dp = opendir(dir);

  if (dp == NULL) {
    printf("Unable to open %s directory\n", dir_p);
    return -1;
  }

  while ((ep = readdir(dp)) != NULL) {
    if (strncmp(base_p, ep->d_name, strlen(base_p)) == 0) {
      int i = strtol(ep->d_name + strlen(base_p), NULL, 10);

      if ((*first_bin < 0) || (i < *first_bin)) {
        *first_bin = i;
      }
      if ((*last_bin < 0) || (i > *last_bin)) {
        *last_bin = i;
      }
    }
  }

  closedir(dp);

  return 0;
}

//
// If estimated bin == initial_bin, return
// If estimated bin == previous bin, return higher bin (in between)

int get_bin_from_time(char *prefix, int initial_bin, int previous_bin,
                      int first_bin, int last_bin, time_t target,
                      int *number_subints_per_file) {
  int status;
  char filename[128], header[256], nfd[32];
  FILE *file;
  int nch;
  double freq;
  double samp_rate;
  float length;
  struct tm mjd_tm = {0};
  time_t mjd;

  sprintf(filename, "%s_%06i.bin", prefix, initial_bin);

  file = fopen(filename, "r");

  if (file == NULL) {
    return -1;
  }

  status = fread(header, sizeof(char), 256, file);

  fclose(file);

  status =
      sscanf(header,
             "HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf "
             "Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\n",
             nfd, &freq, &samp_rate, &length, &nch, number_subints_per_file);

  strptime(nfd, "%Y-%m-%dT%T", &mjd_tm);
  mjd = mktime(&mjd_tm);

  int estimated_bin =
      floor(initial_bin +
            difftime(target, mjd) / (length * *number_subints_per_file));

  // Clip to file range
  if (estimated_bin < first_bin) {
    estimated_bin = first_bin;
  }

  if (estimated_bin > last_bin) {
    estimated_bin = last_bin;
  }

  if (estimated_bin == initial_bin) {
    return estimated_bin;
  } else if (estimated_bin == previous_bin) {
    // Prevent infinite ing poin
    // Return higher as the target is in the gap between
    return estimated_bin > initial_bin ? estimated_bin : initial_bin;
  }

  return get_bin_from_time(prefix, estimated_bin, initial_bin, first_bin,
                           last_bin, target, number_subints_per_file);
}

int get_subs_from_datestrings(char *prefix, char *start, char *end,
                              int *start_bin, int *num_integrations) {
  int k, status;
  char filename[128], header[256], nfd[32];
  FILE *file;
  int number_subints_per_file;
  int nch;
  double freq;
  double samp_rate;
  float length;
  int first_bin = -1;
  int last_bin = -1;

  // Get avaiable file range to constrain search
  get_bin_range(prefix, &first_bin, &last_bin);

  if (start != NULL) {
    struct tm tm = {0};
    time_t t;
    strptime(start, "%Y-%m-%dT%T", &tm);
    t = mktime(&tm);

    int s = get_bin_from_time(prefix, first_bin, first_bin, first_bin, last_bin,
                              t, &number_subints_per_file);

    if (s < 0) {
      return -1;
    }

    *start_bin = s;
  }

  if (end != NULL) {
    struct tm tm = {0};
    time_t t;
    strptime(end, "%Y-%m-%dT%T", &tm);
    t = mktime(&tm);

    int s = get_bin_from_time(prefix, first_bin, first_bin, first_bin, last_bin,
                              t, &number_subints_per_file);

    if (s < *start_bin) {
      return -1;
    }

    *num_integrations = (s - *start_bin + 1) * number_subints_per_file;
  }

  return 0;
}

struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin,double foff)
{
  int i,j,k,l,status,msub,ibin,nadd,nbits=-32;
  char filename[128],header[256],nfd[32];
  FILE *file;
  struct spectrogram s;
  float *z,zavg,zstd;
  char *cz;
  int nch,j0,j1;
  double freq,samp_rate;
  float length;
  int nchan,dummy;

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
  if (strstr(header,"NBITS         8")==NULL) {
    status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\n",s.nfd0,&s.freq,&s.samp_rate,&length,&nch,&msub);
  } else {
    status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nNBITS         8\nMEAN         %f\nRMS          %f",s.nfd0,&s.freq,&s.samp_rate,&length,&nch,&msub,&zavg,&zstd);
    nbits=8;
  }
  s.freq+=foff;

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

  // Read whole file if not specified
  if (nsub==0 && msub>0)
    nsub=msub;

  // Number of subints
  s.nsub=nsub/nbin;
  s.msub=msub;
  s.isub=isub;

  printf("Allocating %.2f MB of memory\n",(4* (float) s.nchan * (float) s.nsub)/(1024 * 1024));

  // Allocate
  s.z=(float *) malloc(sizeof(float)*s.nchan*s.nsub);
  s.zavg=(float *) malloc(sizeof(float)*s.nsub);
  s.zstd=(float *) malloc(sizeof(float)*s.nsub);
  z=(float *) malloc(sizeof(float)*nch);
  cz=(char *) malloc(sizeof(char)*nch);
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
      if (nbits==-32)
	status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\n",nfd,&freq,&samp_rate,&length,&nchan,&msub);
      else if (nbits==8)
	status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nNBITS         8\nMEAN         %f\nRMS          %f",nfd,&freq,&samp_rate,&length,&nchan,&dummy,&zavg,&zstd);

      s.mjd[i]+=nfd2mjd(nfd)+0.5*length/86400.0;
      s.length[i]+=length;
      nadd++;

      // Read buffer
      if (nbits==-32) {
	status=fread(z,sizeof(float),nch,file);
      } else if (nbits==8) {
	status=fread(cz,sizeof(char),nch,file);
	for (j=0;j<nch;j++)
	  z[j]=6.0/256.0*(float) cz[j]*zstd+zavg;
      }
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
  if(nadd>0)
  {
    s.mjd[i]/=(float) nadd;

    for (j=0;j<s.nchan;j++)
      s.z[i+s.nsub*j]/=(float) nadd;
  }

  // Swap frequency range
  if (f0>0.0 && df0>0.0) {
    s.freq=f0;
    s.samp_rate=df0;
  }

  // Compute averages
  for (i=0;i<s.nsub;i++) {
    s.zavg[i]=0.0;
    for (j=0;j<s.nchan;j++)
      if (!isnan(s.z[i+s.nsub*j]) && !isinf(s.z[i+s.nsub*j]))
	s.zavg[i]+=s.z[i+s.nsub*j];
    s.zavg[i]/=(float) s.nchan;
  }

  // Compute deviations
  for (i=0;i<s.nsub;i++) {
    s.zstd[i]=0.0;
    for (j=0;j<s.nchan;j++)
      if (!isnan(s.z[i+s.nsub*j]) && !isinf(s.z[i+s.nsub*j]))
	s.zstd[i]+=pow(s.zavg[i]-s.z[i+s.nsub*j],2);
    s.zstd[i]=sqrt(s.zstd[i]/(float) s.nchan);
  }

  // Compute limits
  for (i=0;i<s.nsub;i++) {
    if (i==0) {
      s.zmin=s.zavg[i]-1.0*s.zstd[i];
      s.zmax=s.zavg[i]+1.0*s.zstd[i];
    } else {
      if (s.zavg[i]-1.0*s.zstd[i]<s.zmin) s.zmin=s.zavg[i]-1.0*s.zstd[i];
      if (s.zavg[i]+1.0*s.zstd[i]>s.zmax) s.zmax=s.zavg[i]+1.0*s.zstd[i];
    }
  }

  double z1,z2;
  zscale(&s, s.nsub, 0.25,&z1, &z2);
  printf("z1 = %f, z2 = %f\n", z1, z2);
  s.zmin = z1;
  s.zmax = z2;
  // Free
  free(z);
  free(cz);

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

void free_spectrogram(struct spectrogram s)
{
  free(s.z);
  free(s.zavg);
  free(s.zstd);
  free(s.mjd);
  free(s.length);
}
