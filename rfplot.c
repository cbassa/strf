#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>
#include "rftime.h"
#include "rfio.h"
#include "rftrace.h"

#define LIM 128
#define NMAX 64

struct select {
  int flag,n;
  float x[NMAX],y[NMAX],w;
};

void dec2sex(double x,char *s,int f,int len);
void time_axis(double *mjd,int n,float xmin,float xmax,float ymin,float ymax);
void usage(void);
void plot_traces(struct trace *t,int nsat);
struct trace fit_trace(struct spectrogram s,struct select sel,int site_id,int graves);
void convolve(float *y,int n,float *w,int m,float *z);
float gauss(float x,float w);
void quadfit(float x[],float y[],int n,float a[]);

// Fit trace
struct trace locate_trace(struct spectrogram s,struct trace t,int site_id)
{
  int i,j,k,l,sn,w=100.0;
  int i0,i1,j0,j1,jmax;
  double f,fmin;
  float x,y,s1,s2,z,za,zs,zm,sigma;
  FILE *file;
  char filename[64];
  
  sprintf(filename,"track_%05d_%08.3f.dat",t.satno,t.freq0);

  // Open file
  file=fopen(filename,"a");

  fmin=(s.freq-0.5*s.samp_rate)*1e-6;

  // Loop over trace
  for (i=0;i<t.n;i++) {
    // Skip when satellite is below the horizon
    if (t.za[i]>90.0)
      continue;
    
    // Compute position
    y=(t.freq[i]-fmin)*s.nchan/(s.samp_rate*1e-6);
    j0=(int) floor(y-w);
    j1=(int) floor(y+w);

    // Keep in range
    if (j0<0)
      j0=0;
    if (j1>=s.nchan)
      j1=s.nchan;
    
    // Find maximum and significance
    zm=0.0;
    jmax=0;
    s1=0.0;
    s2=0.0;
    sn=0;
    for (j=j0;j<j1;j++) {
      z=s.z[i+s.nsub*j];
      s1+=z;
      s2+=z*z;
      sn++;
      if (z>zm) {
	zm=z;
	jmax=j;
	}
    }
    za=s1/(float) sn;
    zs=sqrt(s2/(float) sn-za*za);
    sigma=(zm-za)/zs;

    // Store
    if (sigma>5.0 && s.mjd[i]>1.0) {
      f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
      fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sigma,site_id);
      cpgpt1((float) i,(float) jmax,17);
    }
  }

  // Close file
  fclose(file);

  return t;
}


void filter(struct spectrogram s,int site_id)
{
  int i,j,k,l,jmax,zmax;
  float s1,s2,avg,std,dz;
  FILE *file;
  double f;
  int *mask;
  float sigma=5;

  mask=(int *) malloc(sizeof(int)*s.nchan);

  // Open file
  file=fopen("filter.dat","w");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
    if (s.mjd[i]==0.0)
      continue;

    // Set mask
    for (j=0;j<s.nchan;j++)
      mask[j]=1;

    // Iterate to remove outliers
    for (k=0;k<10;k++) {

      // Find average
      for (j=0,s1=s2=0.0;j<s.nchan;j++) {
	if (mask[j]==1) {
	  s1+=s.z[i+s.nsub*j];
	  s2+=1.0;
	}
      }
      avg=s1/s2;
      
      // Find standard deviation
      for (j=0,s1=s2=0.0;j<s.nchan;j++) {
	if (mask[j]==1) {
	  dz=s.z[i+s.nsub*j]-avg;
	  s1+=dz*dz;
	  s2+=1.0;
	}
      }
      std=sqrt(s1/s2);

      // Update mask
      for (j=0,l=0;j<s.nchan;j++) {
	if (fabs(s.z[i+s.nsub*j]-avg)>sigma*std) {
	  mask[j]=0;
	  l++;
	}
      }
    }
    // Reset mask
    for (j=0;j<s.nchan;j++) {
      if (s.z[i+s.nsub*j]-avg>sigma*std) 
	mask[j]=1;
      else
	mask[j]=0;
    }    
    /*
    // Find maximum when points are adjacent
    for (j=0;j<s.nchan-1;j++) {
      if (mask[j]==1 && mask[j+1]==1) {
	if (s.z[i+s.nsub*j]<s.z[i+s.nsub*(j+1)])
	  mask[j]=0;
      }
    }
    for (j=s.nchan-2;j>=0;j--) {
      if (mask[j]==1 && mask[j-1]==1) {
	if (s.z[i+s.nsub*j]<s.z[i+s.nsub*(j-1)])
	  mask[j]=0;
      }
    }
    */
    // Mark points
    for (j=0;j<s.nchan;j++) {
      if (mask[j]==1) {
	f=s.freq-0.5*s.samp_rate+(double) j*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
	cpgpt1((float) i+0.5,(float) j+0.5,17);
      }
    }
  }
  fclose(file);

  free(mask);

  return;
}

void peakfind(struct spectrogram s,int site_id,int i0,int i1,int j0,int j1)
{
  int i,j,k,l,m=21,n;
  float *w,*y,*sy,*a,*b,*c,d[3],dx[3],dw=1.0,x0,c0=-0.000008;
  double f;
  FILE *file;

  n=j1-j0;

  // Allocate
  y=(float *) malloc(sizeof(float)*n);
  sy=(float *) malloc(sizeof(float)*n);
  a=(float *) malloc(sizeof(float)*n); 
  b=(float *) malloc(sizeof(float)*n); 
  c=(float *) malloc(sizeof(float)*n);

  // Make gaussian smoothing filter
  w=(float *) malloc(sizeof(float)*m);
  for (i=0;i<m;i++) 
    w[i]=gauss((float) (i-m/2),dw);

  // Open file
  file=fopen("peakfind.dat","w");

  // Loop over subints
  for (i=i0;i<i1;i++) {
    if (s.mjd[i]==0.0)
      continue;

    // Fill array
    for (j=0;j<n;j++)
      y[j]=s.z[i+s.nsub*(j0+j)];

    // Convolve
    convolve(y,n,w,m,sy);
    
    // Fit parabolas
    dx[0]=-1.0;
    dx[1]=0.0;
    dx[2]=1.0;
    for (j=1;j<n-1;j++) {
      quadfit(dx,&sy[j-1],3,d);
      a[j]=d[0];
      b[j]=d[1];
      c[j]=d[2];
    }

    // Mark points
    for (j=0;j<n-1;j++) {
      if (b[j]>0.0 && b[j+1]<0.0 && c[j]<c0) {
	x0=(float) (j+j0)+b[j]/(b[j]-b[j+1]);
	f=s.freq-0.5*s.samp_rate+(double) x0*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
	cpgpt1((float) i+0.5,x0+0.5,17);
      }
    }
  }

  // Close
  fclose(file);

  // Free
  free(y);
  free(sy);
  free(a);
  free(b);
  free(c);
  free(w);

  return;
}

int main(int argc,char *argv[])
{
  struct spectrogram s;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float cool_l[]={-0.5,0.0,0.17,0.33,0.50,0.67,0.83,1.0,1.7};
  float cool_r[]={0.0,0.0,0.0,0.0,0.6,1.0,1.0,1.0,1.0};
  float cool_g[]={0.0,0.0,0.0,1.0,1.0,1.0,0.6,0.0,1.0};
  float cool_b[]={0.0,0.3,0.8,1.0,0.3,0.0,0.0,0.0,1.0};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float xmin,xmax,ymin,ymax,zmin,zmax=8.0;
  int i,j,k,flag=0,isel=0,sn;
  int redraw=1,mode=0,posn=0,click=0,graves=0,grid=0;
  float dt,zzmax,s1,s2,z,za,sigma,zs,zm;
  int ix=0,iy=0,isub=0;
  int i0,j0,i1,j1,jmax;
  float width=1500;
  float x,y,x0,y0;
  char c;
  char path[128],xlabel[64],ylabel[64],filename[32],tlefile[128];
  int sec,lsec,ssec;
  char stime[16];
  double fmin,fmax,fcen,f;
  FILE *file;
  int arg=0,nsub=3600,nbin=1;
  double f0=0.0,df0=0.0;
  int foverlay=1;
  struct trace *t,tf;
  int nsat,satno,status;
  struct select sel;
  char *env;
  int site_id=0;
  int cmap=0;
  double foff=0.0,mjdgrid=0.0;
  int jj0,jj1;

  // Get site
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    site_id=atoi(env);
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }
  env=getenv("ST_TLEDIR");
  sprintf(tlefile,"%s/bulk.tle",env);  

  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:hc:C:gm:o:"))!=-1) {
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
	
      case 'b':
	nbin=atoi(optarg);
	break;
	
      case 'f':
	f0=(double) atof(optarg);
	break;
	
      case 'w':
	df0=(double) atof(optarg);
	break;
	
      case 'z':
	zmax=atof(optarg);
	break;

      case 'g':
	graves=1;
	break;

      case 'h':
	usage();
	return 0;

      case 'm':
	cmap=atoi(optarg);
	if (cmap>2)
	  cmap=0;
	break;

      case 'o':
	foff=(double) atof(optarg);
	break;

      case 'c':
	strcpy(tlefile,optarg);
	break;
	
      case 'C':
	site_id=atoi(optarg);
	break;

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
  s=read_spectrogram(path,isub,nsub,f0,df0,nbin,foff);
  
  printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

  // Exit on empty data
  if (s.nsub==0)
    return 0;

  // Compute traces
  t=compute_trace(tlefile,s.mjd,s.nsub,site_id,s.freq*1e-6,s.samp_rate*1e-6,&nsat,graves);
  printf("Traces for %d objects for location %d\n",nsat,site_id);

  cpgopen("/xs");
  cpgsch(0.8);
  cpgask(0);

  // Default limits
  xmin=0.0;
  xmax=(float) s.nsub;
  ymin=0.0;
  ymax=(float) s.nchan;
  zmin=0.0;

  // Set trace
  tf.n=0;
  
  // Set selection
  isel=0;
  sel.n=0;
  sel.w=50.0;

  // Forever loop
  for (;;) {
    if (redraw==1) {
      //      cpgeras();
      cpgpage();
      cpgsci(1);
      /*
      cpgsvp(0.1,0.95,0.9,0.95);
      cpgswin(xmin,xmax,s.zmin,s.zmax);
      cpgbox("BCTS",0.,0,"BCTSN",0.,0);

      for (i=0;i<s.nsub;i++) {
	if (i==0)
	  cpgmove((float) i,s.zavg[i]);
	else
	  cpgdraw((float) i,s.zavg[i]);
      }

      cpgsvp(0.1,0.95,0.1,0.85);
      */
      cpgsvp(0.1,0.95,0.1,0.95);
      cpgswin(xmin,xmax,ymin,ymax);

      if (cmap==2) {
	cpggray(s.z,s.nsub,s.nchan,1,s.nsub,1,s.nchan,zmax,zmin,tr);
      } else {
	if (cmap==0)
	  cpgctab(cool_l,cool_r,cool_g,cool_b,9,1.0,0.5);
	else if (cmap==1)
	  cpgctab(heat_l,heat_r,heat_g,heat_b,9,1.0,0.5);
	cpgimag(s.z,s.nsub,s.nchan,1,s.nsub,1,s.nchan,zmin,zmax,tr);
      }

      // Pixel axis
      cpgbox("CTSM1",0.,0,"CTSM1",0.,0);

      // Time axis
      cpgbox("B",0.,0,"",0.,0);
      time_axis(s.mjd,s.nsub,xmin,xmax,ymin,ymax);

      // Freq axis
      fmin=s.freq-0.5*s.samp_rate+ymin*s.samp_rate/(float) s.nchan;
      fmax=s.freq-0.5*s.samp_rate+ymax*s.samp_rate/(float) s.nchan;
      fmin*=1e-6;
      fmax*=1e-6;

      cpgswin(xmin,xmax,fmin,fmax);
      if (foverlay==1) {
	cpgsci(3);
	plot_traces(t,nsat);
	cpgsci(1);
      }

      // Human readable frequency axis
      fcen=0.5*(fmax+fmin);
      fcen=floor(1000*fcen)/1000.0;
      sprintf(ylabel,"Frequency - %.3f MHz",fcen);
      fmin-=fcen;
      fmax-=fcen;
      cpgswin(xmin,xmax,fmin,fmax);
      cpgbox("",0.,0,"BTSN",0.,0);

      sprintf(xlabel,"UT Date: %.10s",s.nfd0);
      cpglab(xlabel,ylabel," ");

      cpgswin(xmin,xmax,ymin,ymax);

      // Plot selection
      if (sel.n>0) {
	cpgsci(7);
	// Plot points
	for (i=0;i<sel.n;i++) 
	  cpgpt1(sel.x[i],sel.y[i],4);
	// Plot upper bound
	for (i=0;i<sel.n;i++) {
	  if (i==0)
	    cpgmove(sel.x[i],sel.y[i]+sel.w);
	  else
	    cpgdraw(sel.x[i],sel.y[i]+sel.w);
	}
	// Plot lower bound
	for (i=0;i<sel.n;i++) {
	  if (i==0)
	    cpgmove(sel.x[i],sel.y[i]-sel.w);
	  else
	    cpgdraw(sel.x[i],sel.y[i]-sel.w);
	}

	cpgsci(1);
      }

      // Plot grid
      if (grid==1) {
	cpgsci(2);
	for (i=0,flag=0;i<s.nsub-1;i++) {
	  dt=86400.0*(s.mjd[i]-mjdgrid);
	  jj1=(int) (floor) (dt/2.4);
	  if (i==0)
	    jj0=jj1;
	  if (jj1-jj0>0.0) {
	    flag=0;
	    jj0=jj1;
	  }
	  if (jj0%2==0)
	    cpgsls(1);
	  else
	    cpgsls(2);
	  if (flag==0) {
	    cpgmove((float) i,ymin);
	    cpgdraw((float) i,ymax);
	    flag=1;
	  }
	}
	cpgsci(1);
	cpgsls(1);
      }
      redraw=0;
    }

    // Get cursor
    cpgband(mode,posn,x0,y0,&x,&y,&c);

    // Quit
    if (c=='q')
      break;

    // Toggle grid
    if (c=='k') {
      if (grid==0)
	grid=1;
      else
	grid=0;
      mjdgrid=s.mjd[(int) floor(x)];
      redraw=1;
    }

    // Track
    if (c=='t') {
      for (i=0;i<nsat;i++) {
	printf("Locating trace for object %05d\n",t[i].satno);
	locate_trace(s,t[i],4171);
      }
    }

    // Select start
    if (c=='s') {
      sel.x[isel]=x;
      sel.y[isel]=y;
      isel++;
      sel.n=isel;
      redraw=1;
      continue;
    }

    if (c=='g')
      filter(s,site_id);

    if (c=='G') {
      i0=(int) floor(xmin);
      i1=(int) ceil(xmax);
      j0=(int) floor(ymin);
      j1=(int) ceil(ymax);
      if (i0<0)
	i0=0;
      if (i1>=s.nsub)
	i1=s.nsub-1;
      if (j0<0)
	j0=0;
      if (j1>=s.nchan)
	j1=s.nchan-1;
      peakfind(s,site_id,i0,i1,j0,j1);
    }

    // Fit
    if (c=='f') {
      tf=fit_trace(s,sel,site_id,graves);
      tf.site=site_id;
      continue;
    }

    // Identify
    if (c=='i') {
      if (graves==0)
	identify_trace(tlefile,tf,0);
      else
	identify_trace_graves(tlefile,tf,0);
      redraw=1;
      continue;
    }

    // Identify
    if (c=='I') {
      printf("Provide satno: ");
      status=scanf("%d",&satno);
      identify_trace(tlefile,tf,satno);
      redraw=1;
      continue;
    }

    // Undo
    if (c=='u') {
      isel--;
      sel.n=isel;
      redraw=1;
    }

    // Increase
    if (c=='v') {
      zmax*=1.01;
      redraw=1;
      printf("Zmax: %g\n",zmax);
      continue;
    }
    if (c=='b') {
      zmax*=0.99;
      redraw=1;
      printf("Zmax: %g\n",zmax);
      continue;
    }
    
    // Locate
    if (c=='l') {
      i0=(int) x;
      jmax=(int) y;
      for (i=i0;;i++) {
	j0=(int) floor(jmax-10);
	j1=(int) ceil(jmax+10);
	if (i<(int) xmin)
	  break;
	if (i>=(int) xmax)
	  break;
	if (j0<0)
	  j0=0;
	if (j1>=s.nchan)
	  j1=s.nchan-1;
	
	zzmax=0.0;
	jmax=0;
	for (j=j0;j<j1;j++) {
	  if (s.z[i+s.nsub*j]>zzmax) {
	    zzmax=s.z[i+s.nsub*j];
	    jmax=j;
	  }
	}
	printf("%d\n",jmax);
	cpgpt1((float) i,(float) jmax,17);
      }
      i0=(int) x;
      jmax=(int) y;
      for (i=i0;;i--) {
	j0=(int) floor(jmax-10);
	j1=(int) ceil(jmax+10);
	if (i<(int) xmin)
	  break;
	if (i>=(int) xmax)
	  break;
	if (j0<0)
	  j0=0;
	if (j1>=s.nchan)
	  j1=s.nchan-1;
	
	zzmax=0.0;
	jmax=0;
	for (j=j0;j<j1;j++) {
	  if (s.z[i+s.nsub*j]>zzmax) {
	    zzmax=s.z[i+s.nsub*j];
	    jmax=j;
	  }
	}
	printf("%d\n",jmax);
	cpgpt1((float) i,(float) jmax,17);
      }
      continue;
    }

    // Mark single point
    if (c=='D') {
      file=fopen("mark.dat","a");
      i=(int) floor(x);
      j=(int) floor(y);
      f=s.freq-0.5*s.samp_rate+(double) j*s.samp_rate/(double) s.nchan;
      if (s.mjd[i]>1.0) {
	if (graves==0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
	else 
	  fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
	printf("%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
      }
      fclose(file);
    }

    // Mark
    if (c=='m') {
      i0=(int) floor(xmin);
      i1=(int) ceil(xmax);
      j0=(int) floor(ymin);
      j1=(int) ceil(ymax);
      if (i0<0)
	i0=0;
      if (i1>=s.nsub)
	i1=s.nsub-1;
      if (j0<0)
	j0=0;
      if (j1>=s.nchan)
	j1=s.nchan-1;

      file=fopen("out.dat","w");
      // Loop over image
      for (i=i0;i<i1;i++) {
	zzmax=0.0;
	jmax=0;
	s1=0.0;
	s2=0.0;
	sn=0;
	for (j=j0;j<j1;j++) {
	  z=s.z[i+s.nsub*j];
	  if (z>zzmax) {
	    zzmax=z;
	    jmax=j;
	  }
	  s1+=z;
	  s2+=z*z;
	  sn++;
	}
	za=s1/(float) sn;
	zs=sqrt(s2/(float) sn-za*za);
	sigma=(zzmax-za)/zs;

	f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
	if (sigma>5.0 && s.mjd[i]>1.0) {
	  if (graves==0)
	    fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,zzmax,site_id);
	  else
	    fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,zzmax,site_id);
	  cpgpt1((float) i,(float) jmax,17);
	}
      }
      fclose(file);
    }

    // Mark
    if (c=='a') {
      i0=(int) floor(xmin);
      i1=(int) ceil(xmax);
      j0=(int) floor(ymin);
      j1=(int) ceil(ymax);
      if (i0<0)
	i0=0;
      if (i1>=s.nsub)
	i1=s.nsub-1;
      if (j0<0)
	j0=0;
      if (j1>=s.nchan)
	j1=s.nchan-1;

      printf("Provide filename: ");
      status=scanf("%s",filename);
      
      file=fopen(filename,"a");
      // Loop over image
      for (i=i0;i<i1;i++) {
	zzmax=0.0;
	jmax=0;
	for (j=j0;j<j1;j++) {
	  if (s.z[i+s.nsub*j]>zzmax) {
	    zzmax=s.z[i+s.nsub*j];
	    jmax=j;
	  }
	}
	f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,zzmax,site_id);
	cpgpt1((float) i,(float) jmax,17);
      }
      fclose(file);
    }

    // Center
    if (c=='c') {
      xmin=x-width;
      xmax=x+width;
      ymin=y-width;
      ymax=y+width;
      redraw=1;
      continue;
    }

    // Color map
    if (c=='C') {
      cmap++;
      if (cmap>2)
	cmap=0;
      redraw=1;
      continue;
    }

    // Toggle overlay
    if (c=='p' || c=='X') {
      if (foverlay==0)
	foverlay=1;
      else if (foverlay==1)
	foverlay=0;
      redraw=1;
    }

    // Width
    if (isdigit(c)) {
      width=1000.0/(c-'0');
      xmin=x-width;
      xmax=x+width;
      ymin=y-width;
      ymax=y+width;
      redraw=1;
      continue;
    }

    // Zoom
    if (c=='+' || c=='=') {
      width/=1.5;
      xmin=x-width;
      xmax=x+width;
      ymin=y-width;
      ymax=y+width;
      redraw=1;
      continue;
    }

    // Unzoom
    if (c=='x' || c=='-') {
      width*=1.5;
      xmin=x-width;
      xmax=x+width;
      ymin=y-width;
      ymax=y+width;
      redraw=1;
      continue;
    }

    // Recompute traces
    if (c=='R') {
      t=compute_trace(tlefile,s.mjd,s.nsub,site_id,s.freq*1e-6,s.samp_rate*1e-6,&nsat,graves);
      redraw=1;
      continue;
    }

    // Reset
    if (c=='r') {
      xmin=0.0;
      xmax=(float) s.nsub;
      ymin=0.0;
      ymax=(float) s.nchan;
      isel=0;
      sel.n=0;
      redraw=1;
      continue;
    }

    // Zoom
    if (c=='z') {
      click=1;
      mode=2;
    }

    // Pan
    if (c=='\t') {

      // Set area
      x=width*(ix+0.5);
      y=width*(iy+0.5);
      xmin=x-0.75*width;
      xmax=x+0.75*width;
      ymin=y-0.75*width;
      ymax=y+0.75*width;

      // Increment
      iy++;
      if (width*ix>(float) s.nsub) {
	ix=0;
	iy=0;
      }
      if (width*iy>(float) s.nchan) {
	ix++;
	iy=0;
      }
      redraw=1;
      continue;
    }

        // Execute zoom, or box delete
    if (c=='A') {
      if (click==0) {
	click=1;
      } else if (click==1 && mode==2) {
	xmin=(x0<x) ? x0 : x;
	xmax=(x0>x) ? x0 : x;
	ymin=(y0<y) ? y0 : y;
	ymax=(y0>y) ? y0 : y;

	click=0;
	mode=0;
	redraw=1;
      } else {
	click=0;
	mode=0;
	redraw=1;
      }
    }

    // Save
    x0=x;
    y0=y;
  }

  cpgend();

  // Free
  free(s.z);
  free(s.zavg);
  free(s.zstd);
  free(s.mjd);
  if (tf.n>0) {
    free(tf.mjd);
    free(tf.freq);
    free(tf.za);
  }
  for (i=0;i<nsat;i++) {
    free(t[i].mjd);
    free(t[i].freq);
    free(t[i].za);
  }

  return 0;
}

// Convert Decimal into Sexagesimal
void dec2sex(double x,char *s,int f,int len)
{
  int i;
  double sec,deg,min;
  char sign;
  char *form[]={"::",",,","hms","  "};

  sign=(x<0 ? '-' : ' ');
  x=3600.*fabs(x);

  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  //  deg=fmod(x,60.);
  deg=x;

  if (len==7) sprintf(s,"%c%02i%c%02i%c%07.4f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==6) sprintf(s,"%c%02i%c%02i%c%06.3f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==5) sprintf(s,"%c%02i%c%02i%c%05.2f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==4) sprintf(s,"%c%02i%c%02i%c%04.1f%c",sign,(int) deg,form[f][0],(int) min,form[f][1],sec,form[f][2]);
  if (len==2) sprintf(s,"%c%02i%c%02i%c%02i%c",sign,(int) deg,form[f][0],(int) min,form[f][1],(int) floor(sec),form[f][2]);

  return;
}

void time_axis(double *mjd,int n,float xmin,float xmax,float ymin,float ymax)
{
  int i,imin,imax;
  double mjdt,mjdmin,mjdmax;
  float dt,t,tmin,tmax;
  int lsec,ssec,sec;
  char stime[16];

  // Find extrema
  for (i=0;i<n;i++) {
    if (i==0) {
      mjdmin=mjd[i];
      mjdmax=mjd[i];
    } else {
      if (mjd[i]>mjdmax) mjdmax=mjd[i];
    }
  }
  dt=(float) 86400*(mjdmax-mjdmin);

  // Choose tickmarks
  if (dt>43000) {
    lsec=10800;
    ssec=3600;
  } else if (dt>21600) {
    lsec=10800;
    ssec=3600;
  } else if (dt>7200) {
    lsec=1800;
    ssec=300;
  } else if (dt>3600) {
    lsec=600;
    ssec=120;
  } else if (dt>900) {
    lsec=300;
    ssec=60;
  } else {
    lsec=60;
    ssec=10;
  }

  // Extrema
  tmin=86400.0*(mjdmin-floor(mjdmin));
  tmax=tmin+dt;
  tmin=lsec*floor(tmin/lsec);
  tmax=lsec*ceil(tmax/lsec);

  // Large tickmarks
  for (t=tmin;t<=tmax;t+=lsec) {
    mjdt=floor(mjdmin)+t/86400.0;
    if (mjdt>=mjdmin && mjdt<mjdmax) {
      for (i=0;i<n-1;i++)
	if (mjdt>=mjd[i] && mjdt<mjd[i+1])
	  break;
      sec=(int) floor(fmod(t,86400.0));
      dec2sex(((float) sec+0.1)/3600.0,stime,0,2);
      stime[6]='\0';
      cpgtick(xmin,ymin,xmax,ymin,((float) i-xmin)/(xmax-xmin),0.5,0.5,0.3,0.0,stime);
    }
  }

  // Small tickmarks
  for (t=tmin;t<=tmax;t+=ssec) {
    mjdt=floor(mjdmin)+t/86400.0;
    if (mjdt>=mjdmin && mjdt<mjdmax) {
      for (i=0;i<n-1;i++)
	if (mjdt>=mjd[i] && mjdt<mjd[i+1])
	  break;
      sec=(int) floor(t);
      cpgtick(xmin,ymin,xmax,ymin,((float) i-xmin)/(xmax-xmin),0.25,0.25,1.0,1.0,"");
    }
  }

  return;
}

void usage(void)
{
  printf("rfplot: plot RF observations\n\n");
  printf("-p <path>    Input path to file /a/b/c_??????.bin\n");
  printf("-s <start>   Number of starting subintegration [0]\n");
  printf("-l <length>  Number of subintegrations to plot [3600]\n");
  printf("-b <nbin>    Number of subintegrations to bin [1]\n");
  printf("-z <zmax>    Image scaling upper limit [8.0]\n");
  printf("-f <freq>    Frequency to zoom into (Hz)\n");
  printf("-w <bw>      Bandwidth to zoom into (Hz)\n");
  printf("-C <site>    Site ID\n");
  printf("-c <catalog> TLE catalog\n");
  printf("-g           GRAVES data\n");
  printf("-h           This help\n");

  return;
}

void plot_traces(struct trace *t,int nsat)
{
  int i,j,flag,textflag;
  char text[8];

  // Loop over objects
  for (i=0;i<nsat;i++) {
    sprintf(text," %d",t[i].satno);

    // Plot label at start of trace
    if (t[i].za[0]<=90.0)
	cpgtext(0.0,(float) t[i].freq[0],text);

    // Loop over trace
    for (j=0,flag=0,textflag=0;j<t[i].n;j++) {
      // Plot label for rising sources
      if (j>0 && t[i].za[j-1]>90.0 && t[i].za[j]<=90.0)
	cpgtext((float) j,(float) t[i].freq[j],text);

      // Plot line
      if (flag==0) {
	cpgmove((float) j,t[i].freq[j]);
	flag=1;
      } else {
	cpgdraw((float) j,t[i].freq[j]);
      }

      // Below horizon
      if (t[i].za[j]>90.0)
	flag=0;
      else
	flag=1;
    }	  
  }

  return;
}

// Fit trace
struct trace fit_trace(struct spectrogram s,struct select sel,int site_id,int graves)
{
  int i,j,k,l,sn;
  int i0,i1,j0,j1,jmax;
  double f;
  float x,y,s1,s2,z,za,zs,zm,sigma;
  struct trace t;
  FILE *file;

  // Set up trace
  t.satno=99999;
  t.n=(int) ceil(sel.x[sel.n-1]-sel.x[0]);
  t.mjd=(double *) malloc(sizeof(double)*t.n);
  t.freq=(double *) malloc(sizeof(double)*t.n);
  t.za=(float *) malloc(sizeof(float)*t.n);

  // Open file
  file=fopen("out.dat","w");

  // Loop over selected regions
  for (k=0,l=0;k<sel.n-1;k++) {
    for (x=sel.x[k];x<=sel.x[k+1];x+=1.0) {
      y=(x-sel.x[k])/(sel.x[k+1]-sel.x[k])*(sel.y[k+1]-sel.y[k])+sel.y[k];
      i=(int) floor(x);
      j0=(int) floor(y-sel.w);
      j1=(int) floor(y+sel.w);

      // Keep in range
      if (j0<0)
	j0=0;
      if (j1>=s.nchan)
	j1=s.nchan;

      // Find maximum and significance
      zm=0.0;
      jmax=0;
      s1=0.0;
      s2=0.0;
      sn=0;
      for (j=j0;j<j1;j++) {
	z=s.z[i+s.nsub*j];
	s1+=z;
	s2+=z*z;
	sn++;
	if (z>zm) {
	  zm=z;
	  jmax=j;
	}
      }
      za=s1/(float) sn;
      zs=sqrt(s2/(float) sn-za*za);
      sigma=(zm-za)/zs;

      // Store
      if (sigma>5.0 && s.mjd[i]>1.0) {
	f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
	if (graves==0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sigma,site_id);
	else
	  fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,sigma,site_id);
	cpgpt1((float) i,(float) jmax,17);
	t.mjd[l]=s.mjd[i];
	t.freq[l]=f;
	t.za[l]=0.0;
	l++;
      }
    }
  }
  t.n=l;

  // Close file
  fclose(file);

  return t;
}

// Discrete convolution
void convolve(float *y,int n,float *w,int m,float *z)
{
  int i,j,k,imid;

  imid=m/2;
  for (i=0;i<n;i++) {
    z[i]=0.0;
    for (j=0;j<m;j++) {
      k=i-imid+j;
      if (k<0 || k>n-1)
	continue;
      z[i]+=w[j]*y[k];
    }
  }

  return;
}

// Gaussian
float gauss(float x,float w)
{
  float c;

  c=pow(x/w,2);

  return exp(-0.5*c)/(sqrt(2.0*M_PI)*w);
}

// Quadratic fit
void quadfit(float x[],float y[],int n,float a[])
{
  int i;
  float p,q,r,s,t,u,v,d;

  p=q=r=s=t=u=v=0.;
  for (i=0;i<n;i++) {
    p+=x[i];
    q+=pow(x[i],2);
    r+=pow(x[i],3);
    s+=pow(x[i],4);
    t+=y[i];
    u+=x[i]*y[i];
    v+=pow(x[i],2)*y[i];
  }

  d=n*q*s+2.*p*q*r-q*q*q-p*p*s-(float) n*r*r;

  a[0]=(q*s*t+q*r*u+p*r*v-q*q*v-p*s*u-r*r*t)/d;
  a[1]=((float) n*s*u+p*q*v+q*r*t-q*q*u-p*s*t-(float) n*r*v)/d;
  a[2]=((float) n*q*v+p*r*t+p*q*u-q*q*t-p*p*v-(float) n*r*u)/d;

  return;
}
