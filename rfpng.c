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
struct trace locate_trace(struct spectrogram s,struct select sel,int site_id);

int main(int argc,char *argv[])
{
  struct spectrogram s;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float cool_l[]={-0.5,0.0,0.17,0.33,0.50,0.67,0.83,1.0,1.7};
  float cool_r[]={0.0,0.0,0.0,0.0,0.6,1.0,1.0,1.0,1.0};
  float cool_g[]={0.0,0.0,0.0,1.0,1.0,1.0,0.6,0.0,1.0};
  float cool_b[]={0.0,0.3,0.8,1.0,0.3,0.0,0.0,0.0,1.0};
  float xmin,xmax,ymin,ymax,zmin,zmax=8.0;
  int i,j,k;
  float dt,zzmax,s1,s2;
  int ix=0,iy=0,isub=0;
  int i0,j0,i1,j1,jmax;
  float width=1500;
  float x,y,x0,y0;
  char c;
  char path[128],xlabel[64],ylabel[64],filename[32],tlefile[128],pngfile[128];
  int sec,lsec,ssec;
  char stime[16];
  double fmin,fmax,fcen,f;
  FILE *file;
  int arg=0,nsub=900,nbin=1;
  double f0=0.0,df0=0.0,dy=2500;
  int foverlay=1;
  struct trace *t,tf;
  int nsat,satno;
  struct select sel;
  char *env;
  int site_id=0;

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
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:hc:C:"))!=-1) {
      switch (arg) {
	
      case 'p':
	strcpy(path,optarg);
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

      case 'c':
	strcpy(tlefile,optarg);
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

  for (isub=0;;isub+=15) {
    // Read data
    s=read_spectrogram(path,isub,nsub,f0,df0,nbin);
    if (s.mjd[0]<54000)
      break;

    // Output filename
    sprintf(pngfile,"%.19s_%8.3f.png/png",s.nfd0,s.freq*1e-6);

    printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

    // Compute traces
    t=compute_trace(tlefile,s.mjd,s.nsub,site_id,s.freq*1e-6,s.samp_rate*1e-6,&nsat);
    printf("Traces for %d objects for location %d\n",nsat,site_id);

    printf("%s\n",pngfile);
    cpgopen(pngfile);
    cpgctab(cool_l,cool_r,cool_g,cool_b,9,1.0,0.5);
    cpgsch(0.8);
    cpgask(1);
    
    // Default limits
    xmin=0.0;
    xmax=(float) s.nsub;
    ymin=0;
    ymax=(float) s.nchan;
    zmin=0.0;
    
    cpgsci(1);
    cpgpage();
    cpgsvp(0.1,0.95,0.1,0.95);
    cpgswin(xmin,xmax,ymin,ymax);
  
    cpgimag(s.z,s.nsub,s.nchan,1,s.nsub,1,s.nchan,zmin,zmax,tr);
    
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
    
    // Plot traces
    cpgswin(xmin,xmax,fmin,fmax);
    cpgsci(3);
    plot_traces(t,nsat);
    cpgsci(1);
    
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
   
    cpgend();

    // Free
    free(s.z);
    free(s.mjd);

    for (i=0;i<nsat;i++) {
      free(t[i].mjd);
      free(t[i].freq);
      free(t[i].za);
    }
    //    free(tf.mjd);
    //    free(tf.freq);
    //    free(tf.za);

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
  // Override
  lsec=60;
  ssec=10;

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
	cpgmove((float) j,(float) t[i].freq[j]);
	flag=1;
      } else {
	cpgdraw((float) j,(float) t[i].freq[j]);
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

// Locate trace
struct trace locate_trace(struct spectrogram s,struct select sel,int site_id)
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
	fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sigma,site_id);
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
