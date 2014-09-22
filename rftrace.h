struct trace {
  int satno,n,site;
  double *mjd;
  double *freq;
  float *za;
};
struct trace *compute_trace(char *tlefile,double *mjd,int n,int site_id,float fmin,float fmax,int *nsat);
void identify_trace(char *tlefile,struct trace t,int satno);
