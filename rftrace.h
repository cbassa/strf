struct trace {
  char satname[25];
  int satno,n,site,classfd,graves;
  double *mjd;
  double *freq,freq0;
  float *za;
};
struct trace *compute_trace(char *tlefile,double *mjd,int n,int site_id,float fmin,float fmax,int *nsat,int graves,char *freqlist);
void identify_trace(char *tlefile,struct trace t,int satno,char *freqlist);
void identify_trace_graves(char *tlefile,struct trace t,int satno,char *freqlist);
void compute_doppler(char *tlefile,double *mjd,int n,int site_id,int satno,int graves, int skiphigh,char *outfname);
int fgetline(FILE *file,char *s,int lim);
