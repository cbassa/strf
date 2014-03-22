#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>

#define LIM 128
struct spectrogram {
  int nsub,nchan;
  double *mjd;
  double freq,samp_rate;
  float *length;
  float *z;
  char nfd0[32];
};
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);
void dec2sex(double x,char *s,int f,int len);
struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin);
void time_axis(double *mjd,int n,float xmin,float xmax,float ymin,float ymax);
void usage(void);

int main(int argc,char *argv[])
{
  struct spectrogram s;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float cool_l[]={-0.5,0.0,0.17,0.33,0.50,0.67,0.83,1.0,1.7};
  float cool_r[]={0.0,0.0,0.0,0.0,0.6,1.0,1.0,1.0,1.0};
  float cool_g[]={0.0,0.0,0.0,1.0,1.0,1.0,0.6,0.0,1.0};
  float cool_b[]={0.0,0.3,0.8,1.0,0.3,0.0,0.0,0.0,1.0};
  float xmin,xmax,ymin,ymax,zmin,zmax=8.0;
  int i,j,k,flag=0;
  int redraw=1,mode=0,posn=0,click=0;
  float dt,zzmax,s1,s2;
  int ix=0,iy=0,isub=0;
  int i0,j0,i1,j1,jmax;
  float width=500;
  float x,y,x0,y0;
  char c;
  char path[128],xlabel[64],ylabel[64];
  int sec,lsec,ssec;
  char stime[16];
  double fmin,fmax,fcen,f;
  FILE *file;
  int arg=0,nsub=3600,nbin=1;
  double f0=0.0,df0=0.0;

  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:h"))!=-1) {
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

      case 'h':
	usage();

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
  s=read_spectrogram(path,isub,nsub,f0,df0,nbin);

  printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

  cpgopen("/xs");
  cpgctab(cool_l,cool_r,cool_g,cool_b,9,1.0,0.5);
  cpgsch(0.8);

  // Default limits
  xmin=0.0;
  xmax=(float) s.nsub;
  ymin=0.0;
  ymax=(float) s.nchan;
  zmin=0.0;
  
  // Forever loop
  for (;;) {
    if (redraw==1) {
      cpgeras();
      cpgsci(1);

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

      cpgswin(xmin,xmax,fmin,fmax);

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

      redraw=0;
    }

    // Get cursor
    cpgband(mode,posn,x0,y0,&x,&y,&c);

    // Quit
    if (c=='q')
      break;

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
	fprintf(file,"%lf %lf %f 4171\n",s.mjd[i],f,s.z[i+s.nsub*j]);
	printf("%lf %lf %f 4171\n",s.mjd[i],f,s.z[i+s.nsub*j]);
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
	for (j=j0;j<j1;j++) {
	  if (s.z[i+s.nsub*j]>zzmax) {
	    zzmax=s.z[i+s.nsub*j];
	    jmax=j;
	  }
	}
	f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0)
	  fprintf(file,"%lf %lf %f 4171\n",s.mjd[i],f,zzmax);
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

    // Reset
    if (c=='r') {
      xmin=0.0;
      xmax=(float) s.nsub;
      ymin=0.0;
      ymax=(float) s.nchan;
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
      xmin=x-width;
      xmax=x+width;
      ymin=y-width;
      ymax=y+width;

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
  free(s.mjd);



  return 0;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min;
  double mjd,dday;
  float sec;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%f",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Compute Julian Day from Date
double date2mjd(int year,int month,double day)
{
  int a,b;
  double jd;

  if (month<3) {
    year--;
    month+=12;
  }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1852 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
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

struct spectrogram read_spectrogram(char *prefix,int isub,int nsub,double f0,double df0,int nbin)
{
  int i,j,k,l,flag=0,status,msub;
  char filename[128],header[256],nfd[32];
  FILE *file;
  struct spectrogram s;
  float *z;
  int nch,j0,j1;
  double freq,samp_rate;
  float length;
  int nchan;

  // Open first file to get number of channels
  sprintf(filename,"%s_%06d.bin",prefix,isub);
	
  // Open file
  file=fopen(filename,"r");
  if (file==NULL) {
    printf("%s does not exist\n",filename);
    return s;
  }

  // Read header
  status=fread(header,sizeof(char),256,file);
  status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\n",s.nfd0,&s.freq,&s.samp_rate,&length,&nch);

  // Close file
  fclose(file);

  // Compute plotting channel
  if (f0>0.0 && df0>0.0) {
    s.nchan=(int) (df0/s.samp_rate*(float) nch);
    
    j0=(int) ((f0-0.5*df0-s.freq+0.5*s.samp_rate)*(float) nch/s.samp_rate);
    j1=(int) ((f0+0.5*df0-s.freq+0.5*s.samp_rate)*(float) nch/s.samp_rate);
    
    if (j0<0 || j1>nch) 
      fprintf(stderr,"Requested frequency range out of limits\n");
  } else {
    s.nchan=nch;
    j0=0;
    j1=s.nchan;
  }

  // Number of subints
  s.nsub=nsub/nbin;

  // Allocate
  s.z=(float *) malloc(sizeof(float)*s.nchan*s.nsub);
  z=(float *) malloc(sizeof(float)*nch);
  s.mjd=(double *) malloc(sizeof(double)*s.nsub);
  s.length=(float *) malloc(sizeof(float)*s.nsub);

  // Initialize
  for (j=0;j<s.nchan*s.nsub;j++)
    s.z[j]=0.0;
  for (j=0;j<s.nsub;j++)
    s.mjd[j]=0.0;

  // Loop over files
  for (k=0,i=0,l=0;l<nsub;k++) {
    // Generate filename
    sprintf(filename,"%s_%06d.bin",prefix,k+isub);

    // Open file
    file=fopen(filename,"r");
    if (file==NULL) {
      printf("%s does not exist\n",filename);
	  break;
    }
    printf("opened %s\n",filename);

    // Loop over contents of file
    for (;l<nsub;l++) {
      // Read header
      status=fread(header,sizeof(char),256,file);
      if (status==0)
	break;
      status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\n",nfd,&freq,&samp_rate,&length,&nchan);
      s.mjd[i]+=nfd2mjd(nfd)+0.5*length/86400.0;
      s.length[i]+=length;

      // Read buffer
      status=fread(z,sizeof(float),nch,file);
      if (status==0)
	break;
      
      // Copy
      for (j=0;j<s.nchan;j++) 
	s.z[i+s.nsub*j]+=z[j+j0];

      // Increment
      if (l%nbin==nbin-1) {
	// Scale
	s.mjd[i]/=(float) nbin;
	for (j=0;j<s.nchan;j++) 
	  s.z[i+s.nsub*j]/=(float) nbin;

	i++;
      }
    }

    // Close file
    fclose(file);
  }

  // Swap frequency range
  if (f0>0.0 && df0>0.0) {
    s.freq=f0;
    s.samp_rate=df0;
  }

  // Free 
  free(z);

  return s;
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
  printf("-h           This help\n");

  return;
}
