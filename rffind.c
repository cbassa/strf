#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>
#include "rftime.h"
#include "rfio.h"
#include <gsl/gsl_multifit.h>

#define LIM 128
#define NMAX 64


void usage(void)
{
}

float fit_polynomial(float *x,float *y,int n,int m,float *a)
{

  int i,j;
  double chisq;
  gsl_matrix *X,*cov;
  gsl_vector *yy,*w,*c;


  X=gsl_matrix_alloc(n,m);
  yy=gsl_vector_alloc(n);
  w=gsl_vector_alloc(n);

  c=gsl_vector_alloc(m);
  cov=gsl_matrix_alloc(m,m);

  // Fill matrices
  for(i=0;i<n;i++) {
    for (j=0;j<m;j++)
      gsl_matrix_set(X,i,j,pow(x[i],j));
    
    gsl_vector_set(yy,i,y[i]);
    gsl_vector_set(w,i,1.0);
  }

  // Do fit
  gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(n,m);
  gsl_multifit_wlinear(X,w,yy,c,cov,&chisq,work);
  gsl_multifit_linear_free(work);

  // Save parameters
  for (i=0;i<m;i++)
    a[i]=gsl_vector_get(c,(i));

  gsl_matrix_free(X);
  gsl_vector_free(yy);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return chisq;
}

void filter(float *x,float *y,int n,int m,float sigma,int *mask)
{
  int i,j,k,nn;
  float *xx,*yy,*a,chi2,ym;
  float rms;

  // Allocate
  a=(float *) malloc(sizeof(float)*m);
  xx=(float *) malloc(sizeof(float)*n);
  yy=(float *) malloc(sizeof(float)*n);

  // Set intial mask
  for (i=0;i<n;i++) {
    if (y[i]<2)
      mask[i]=1;
    else
      mask[i]=0;
  }

  // Iterations
  for (k=0;k<10;k++) {
    // Apply mask
    for (i=0,j=0;i<n;i+=100) {
      if (mask[i]!=1)
	continue;
      xx[j]=x[i];
      yy[j]=y[i];
      j++;
    }
    nn=j;
    
    // Fit polynomial
    chi2=fit_polynomial(xx,yy,nn,m,a);
    
    // Scale
    for (i=0,rms=0.0,nn=0;i<n;i++) {
      for (j=0,ym=0.0;j<m;j++)
	ym+=a[j]*pow(x[i],j);
      yy[i]=y[i]/ym-1.0;
      if (mask[i]==1) {
	rms+=yy[i]*yy[i];
	nn++;
      }
    }
    rms=sqrt(rms/(float) (nn-1));

    // Update mask
    for (i=0;i<n;i++) {
      if (fabs(yy[i])>sigma*rms)
	mask[i]=0;
    }
  }

  // Recompute final mask
  for (i=0;i<n;i++) {
    if (yy[i]>sigma*rms)
      mask[i]=0;
    else
      mask[i]=1;
  }

  // Free
  free(a);
  free(xx);
  free(yy);

  return;
}

int main(int argc,char *argv[])
{
  int i,j,k,l,j0,j1,m=2,n;
  struct spectrogram s;
  char path[128];
  int isub=0,nsub=60;
  char *env;
  int site_id=0;
  float avg,std;
  int arg=0;
  float sigma=5.0;
  float *x,*y;
  int *mask;
  FILE *file;
  double f,f0,df0;
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
    while ((arg=getopt(argc,argv,"p:f:w:s:l:hc:o:"))!=-1) {
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
	
      case 'o':
	strcpy(filename,optarg);
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
  s=read_spectrogram(path,isub,nsub,f0,df0,1);

  x=(float *) malloc(sizeof(float)*s.nchan);
  y=(float *) malloc(sizeof(float)*s.nchan);
  mask=(int *) malloc(sizeof(int)*s.nchan);

  printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

  file=fopen(filename,"w");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
    // Fill array
    for (j=0;j<s.nchan;j++) {
      x[j]=-0.5*s.samp_rate*1e-6+s.samp_rate*1e-6*(float) j/(float) (s.nchan-1);
      y[j]=s.z[i+s.nsub*j];
    }

    // Filter data
    filter(x,y,s.nchan,m,sigma,mask);

    // Output
    for (j=0;j<s.nchan;j++) {
      if (mask[j]==0) {
	f=s.freq-0.5*s.samp_rate+(double) j*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
      }
    }
  }
  fclose(file);

  // Free
  free(s.z);
  free(s.mjd);
  free(x);
  free(y);
  free(mask);

  return 0;
}

