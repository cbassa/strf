#ifndef RFIO_H
#define RFIO_H
struct spectrogram {
  int nsub,nchan,msub,isub;
  double *mjd;
  double freq,samp_rate;
  float *length;
  float *z,*zavg,*zstd;
  float zmin,zmax;
  char nfd0[32];
};
struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin,double foff);
void write_spectrogram(struct spectrogram s,char *prefix);
void free_spectrogram(struct spectrogram s);
#endif
