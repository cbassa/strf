#ifndef RFIO_H
#define RFIO_H

#include <time.h>

struct spectrogram {
  int nsub,nchan,msub,isub;
  double *mjd;
  double freq,samp_rate;
  float *length;
  float *z,*zavg,*zstd;
  float zmin,zmax;
  char nfd0[32];
};

int get_bin_range(char *prefix, int *first_bin, int *last_bin);
int get_bin_from_time(char *prefix, int initial_bin, int previous_bin,
                      int first_bin, int last_bin, time_t target,
                      int *number_subints_per_file);
int get_subs_from_datestrings(char *prefix, char *start, char *end,
                              int *start_bin, int *num_integrations);
struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin,double foff);
void write_spectrogram(struct spectrogram s,char *prefix);
void free_spectrogram(struct spectrogram s);

#endif // RFIO_H
