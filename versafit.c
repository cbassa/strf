// Versatile Fitting Routine
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int OUTPUT=1; // Print output on screen (1 = yes; 0 = no)
int ERRCOMP=0; // Set reduced Chi-Squared to unity (1 = yes; 0 = no)

int dsmin(double **,double *,int,double,double (*func)(double *));
double **simplex(int,double *,double *);
double parabolic_root(double,double,double,double);

// Versafit fitting routine
//
// Inputs:
//    m:      number of datapoints
//    n:      number of parameters
//    a:      parameters
//    da:     expected spread in parameters
//    func:   function to fit (Chi-squared function)
//    dchisq  difference in Chi-squared
//    tol:    tolerance
//    opt:    options
//            - n: no output
void versafit(int m,int n,double *a,double *da,double (*func)(double *),double dchisq,double tol,char *opt)
{
  int i,j,k,l,nfunk,kmax=50;
  double chisqmin;
  double *b,*db;
  double **p,*y;
  double d[2],errcomp;

  // Decode options
  if (strchr(opt,'n')!=NULL) OUTPUT=0;
  if (strchr(opt,'e')!=NULL) ERRCOMP=1;

  // Intialize y
  y=(double *) malloc(sizeof(double) * (n+1));

  if (dchisq>=0.) {
    // Compute simplex and minimize function
    p=simplex(n,a,da);
    nfunk=dsmin(p,y,n,tol,func);

    // Average parameters
    for (i=0;i<n;i++) {
      a[i]=0.;
      for (j=0;j<=n;j++) 
	a[i]+=p[j][i];
      a[i]/=(double) (n+1);
    }

    // Compute minimum
    chisqmin=func(a);

    // Compute error compensation
    if (ERRCOMP) errcomp=sqrt(chisqmin/(double) (m-n));
  }

  // Basic Information
  if (OUTPUT) {
    printf("VersaFIT:\n");
    if (m!=0)
      printf("Number of datapoints: %i\n",m);
    printf("Number of parameters: %i\n",n);
    printf("Chi-squared: %14.5f\n",chisqmin);
    if (m!=0)
      printf("Reduced Chi-squared: %14.5f\n",chisqmin/(double) (m-n));
    if (ERRCOMP) printf("Error compensation: %.4f\n",errcomp);
    printf("Number of iterations: %i\n",nfunk);
  
    printf("\nParameters:\n");

    // No error estimation
    if (dchisq==0.) {
      for (i=0;i<n;i++) 
	printf(" a(%i): %12.5f\n",i+1,a[i]);
    }
  }

  // With error estimation
  if (dchisq!=0.) {
    b=(double *) malloc(sizeof(double) * n);
    db=(double *) malloc(sizeof(double) * n);

    for (i=0;i<n;i++) {
      if (da[i]!=0.) {
	for (j=0;j<n;j++) {
	  b[j]=a[j];
	  db[j]=da[j];
	}
	d[0]=-da[i];
	db[i]=0.;

	for (k=0;k<kmax;k++) {
	  b[i]=a[i]+d[0];

	  // Minimize
	  p=simplex(n,b,db);
	  nfunk+=dsmin(p,y,n,tol,func);

	  // Average parameters
	  for (l=0;l<n;l++) {
	    b[l]=0.;
	    for (j=0;j<=n;j++) 
	      b[l]+=p[j][l];
	    b[l]/=(double) (n+1);
	  }
	  d[0]=parabolic_root(d[0],func(b),chisqmin,dchisq);

	  if (fabs(chisqmin+dchisq-func(b))<tol) break;
	}

	d[1]=-d[0];
	db[i]=0.;

	for (k=0;k<kmax;k++) {
	  b[i]=a[i]+d[1];

	  // Minimize
	  p=simplex(n,b,db);
	  nfunk+=dsmin(p,y,n,tol,func);

	  // Average parameters
	  for (l=0;l<n;l++) {
	    b[l]=0.;
	    for (j=0;j<=n;j++) 
	      b[l]+=p[j][l];
	    b[l]/=(double) (n+1);
	  }
	  d[1]=parabolic_root(d[1],func(b),chisqmin,dchisq);

	  if (fabs(chisqmin+dchisq-func(b))<tol) break;
	}
	da[i]=0.5*(fabs(d[0])+fabs(d[1]));
	if (ERRCOMP) da[i]*=errcomp;
      }
    }

    if (OUTPUT) 
      for (i=0;i<n;i++) 
	printf(" a(%i): %12.5f +- %9.5f\n",i+1,a[i],da[i]);
  }
  if (OUTPUT) printf("\nTotal number of iterations: %i\n",nfunk);

  //  free(p);
  //  free(y);
  //  free(b);
  //  free(db);

  return;
}

// Compute root
double parabolic_root(double x,double y,double y0,double dy)
{
  double a;

  if (fabs(x)<1e-9) {
    printf("Division by zero in function 'parabolic_root'\n");
    x=1e-9;
  }

  a=(y-y0)/(x*x);
  
  return sqrt(fabs(dy/a))*x/fabs(x);
}
