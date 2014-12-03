// Creates Simplex
#include <stdio.h>
#include <stdlib.h>

double **simplex(int n,double *a,double *da)
{
  int i,j;
  double **p;

  // Allocate pointers to rows
  p=(double **) malloc(sizeof(double *) * (n+1));

  // Allocate rows and set pointers
  for (i=0;i<=n;i++)
    p[i]=(double *) malloc(sizeof(double) * (n+1)*n);

  // Fill simplex
  for (i=0;i<=n;i++) {
    for (j=0;j<n;j++) {
      if (i<j) p[i][j]=a[j];
      if (i==j) p[i][j]=a[j]+da[j];
      if (i>j) p[i][j]=a[j]-da[j];
    }
  }  

  return p;
}

