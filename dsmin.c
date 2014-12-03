#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TINY 1.0e-10
#define NMAX 100000
#define SWAP(a,b) {(a)+=(b);(b)=(a)-(b);(a)-=(b);}

// Downhill Simplex Minimization
int dsmin(double **p,double *y,int n,double ftol,double (*func)(double *))
{
  int i,j,nfunk=0;
  int ihi,ilo,ise;
  double *ptry,*pmid,*psum;
  double tol,ytry,rtol,ysave;
  double *vector_sum(double **,int);
  double dsmod(double **,double *,double *,int,double (*func)(double *),int,double);

  // Allocate memory
  psum=(double *) malloc(sizeof(double) * n);

  // Get function values
  for (i=0;i<=n;i++) 
    y[i]=func(p[i]);

  // Sum vectors
  psum=vector_sum(p,n);

  // Start forever loop
  for (;;) {
    // Find high and low point
    ilo=0;
    ihi = (y[0]>y[1]) ? (ise=1,0) : (ise=0,1);
    for (i=0;i<=n;i++) {
      if (y[i]<=y[ilo]) ilo=i;
      if (y[i]>y[ihi]) {
      ise=ihi;
      ihi=i;
      } else if (y[i]>y[ise] && i!=ihi) ise=i;
    }

    // Compute fractional range from highest to lowest point
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

    // Return if fractional tolerance is acceptable
    if (rtol<ftol) 
      break;

    if (nfunk>=NMAX) {
      printf("dsmin: NMAX exceeded!\n");
      return -1;
    }
    nfunk+=2;

    // Reflect simplex
    ytry=dsmod(p,y,psum,n,func,ihi,-1.0);

    if (ytry<=y[ilo]) // Goes right direction, extrapolate by factor 2
      ytry=dsmod(p,y,psum,n,func,ihi,2.0);
    else if (ytry>=y[ise]) { // 1D contraction
      ysave=y[ihi];
      ytry=dsmod(p,y,psum,n,func,ihi,0.5);
      if (ytry>=ysave) {
	for (i=0;i<=n;i++) {
	  if (i!=ilo) {
	    for (j=0;j<n;j++) 
	      p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	    y[i]=(*func)(psum);
	  }
	}
	nfunk+=n;
	
        psum=vector_sum(p,n);
      }
    } else --nfunk;
  }
  free(psum);

  return nfunk;
}

// Sum vectors
double *vector_sum(double **p,int n)
{
  int i,j;
  double sum,*psum;

  psum=(double *) malloc(sizeof(double) * n);

  for (i=0;i<n;i++) {
    sum=0.;
    for (j=0;j<=n;j++)
      sum+=p[j][i];
    psum[i]=sum;
  }

  return psum;
}

// Simplex modification
double dsmod(double **p,double *y,double *psum,int n,double (*func)(double *),int ihi,double fac)
{
  int i;
  double fac1,fac2,ytry,*ptry;

  ptry=(double *) malloc(sizeof(double) * n);
  fac1=(1.0-fac)/(double) n;
  fac2=fac1-fac;

  for (i=0;i<n;i++)
    ptry[i]=psum[i]*fac1-p[ihi][i]*fac2;
  ytry=(*func)(ptry);
  if (ytry<y[ihi]) {
    y[ihi]=ytry;
    for (i=0;i<n;i++) {
      psum[i] += ptry[i]-p[ihi][i];
      p[ihi][i]=ptry[i];
    }
  }

  free(ptry);
  
  return ytry;
}
