struct spectrogram {
  int nsub,nchan;
  double *mjd;
  double freq,samp_rate;
  float *length;
  float *z;
  char nfd0[32];
};
struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin,double foff);
void write_spectrogram(struct spectrogram s,char *prefix);
