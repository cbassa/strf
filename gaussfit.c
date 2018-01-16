#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

struct data {
  size_t n;
  double *x,*y,*sy,*w;
};

int gauss_f(const gsl_vector *q,void *d,gsl_vector *f)
{
  int i,n=((struct data *) d)->n;
  double *x=((struct data *) d)->x;
  double *y=((struct data *) d)->y;
  double a[5],ym,arg,ex,fac;

  // Extract parameters
  for (i=0;i<5;i++)
    a[i]=gsl_vector_get(q,i);

  // Compute values
  for (i=0;i<n;i++) {
    arg=(x[i]-a[0])/a[1];
    ex=exp(-arg*arg);
    fac=2.0*a[2]*ex*arg;
    ym=a[2]*ex+a[3]+a[4]*x[i];
    gsl_vector_set(f,i,ym-y[i]);
  }

  return GSL_SUCCESS;
}

int gauss_df(const gsl_vector *q,void *d,gsl_matrix *J)
{
  int i,n=((struct data *) d)->n;
  double *x=((struct data *) d)->x;
  double *y=((struct data *) d)->y;
  double a[5],ym,arg,ex,fac;

  // Extract parameters
  for (i=0;i<5;i++)
    a[i]=gsl_vector_get(q,i);

  // Compute values
  for (i=0;i<n;i++) {
    arg=(x[i]-a[0])/a[1];
    ex=exp(-arg*arg);
    fac=2.0*a[2]*ex*arg;
    gsl_matrix_set(J,i,0,fac/a[1]);
    gsl_matrix_set(J,i,1,fac*arg/a[1]);
    gsl_matrix_set(J,i,2,ex);
    gsl_matrix_set(J,i,3,1.0);
    gsl_matrix_set(J,i,4,x[i]);
  }

  return GSL_SUCCESS;
}

int fit_gaussian(double *x,double *y,double *sy,int n,double *q,double *sq,int m,double *chisq)
{
  int i;
  const gsl_multifit_nlinear_type *T=gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
  gsl_vector *f;
  gsl_matrix *J;
  struct data d;
  gsl_rng *r;
  int status,info;
  const double xtol=1e-8;
  const double gtol=1e-8;
  const double ftol=0.0;

  // Fill struct
  d.n=n;
  d.x=(double *) malloc(sizeof(double)*d.n);
  d.y=(double *) malloc(sizeof(double)*d.n);
  d.sy=(double *) malloc(sizeof(double)*d.n);
  d.w=(double *) malloc(sizeof(double)*d.n);
  for (i=0;i<d.n;i++) {
    d.x[i]=x[i];
    d.y[i]=y[i];
    d.sy[i]=sy[i];
    d.w[i]=1.0/(d.sy[i]*d.sy[i]);
  }
  
  gsl_vector_view a=gsl_vector_view_array(q,m);
  gsl_vector_view wts=gsl_vector_view_array(d.w,d.n);
  gsl_matrix *covar=gsl_matrix_alloc(m,m);

  // Function definition
  fdf.f=gauss_f;
  fdf.df=gauss_df;
  fdf.fvv=NULL;
  fdf.n=d.n;
  fdf.p=m;
  fdf.params=&d;

  // Allocate and initialize
  w=gsl_multifit_nlinear_alloc(T,&fdf_params,d.n,m);
  gsl_multifit_nlinear_winit(&a.vector,&wts.vector,&fdf,w);
  f=gsl_multifit_nlinear_residual(w);

  // Solve and compute covariance and chi-squared
  status=gsl_multifit_nlinear_driver(20,xtol,gtol,ftol,NULL,NULL,&info,w);
  J=gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar(J,0.0,covar);
  gsl_blas_ddot(f,f,chisq);

  // Extract parameters
  for (i=0;i<m;i++) {
    q[i]=gsl_vector_get(w->x,i);
    sq[i]=sqrt(gsl_matrix_get(covar,i,i));
  }

  // Free
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
  free(d.x);
  free(d.y);
  free(d.sy);
  free(d.w);
  
  return status;
}

