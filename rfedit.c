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
void mjd2nfd(double mjd,char *nfd);

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

void write_spectrogram(struct spectrogram s,char *prefix)
{
  int i,j;
  FILE *file;
  char header[256]="",filename[256],nfd[32];
  float *z;
  double mjd;

  // Allocate
  z=(float *) malloc(sizeof(float)*s.nchan);

  // Generate filename
  sprintf(filename,"%s_%06d.bin",prefix,0);

  // Open file
  file=fopen(filename,"w");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
    // Date
    mjd=s.mjd[i]-0.5*s.length[i]/86400.0;
    mjd2nfd(mjd,nfd);

    // Generate header
    sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nEND\n",nfd,s.freq,s.samp_rate,s.length[i],s.nchan);

    // Copy buffer
    for (j=0;j<s.nchan;j++) 
      z[j]=s.z[i+s.nsub*j];

    // Dump contents
    fwrite(header,sizeof(char),256,file);
    fwrite(z,sizeof(float),s.nchan,file);
  }

  // Close file
  fclose(file);

  // Free
  free(z);

  return;
}

int main(int argc,char *argv[])
{
  struct spectrogram s;
  char prefix[128],outfile[128];
  int arg=0,nsub=3600,nbin=1,isub=0;
  double f0=0.0,df0=0.0;

  // Read arguments
  while ((arg=getopt(argc,argv,"i:o:f:w:s:l:b:z:"))!=-1) {
    switch (arg) {
      
    case 'i':
      strcpy(prefix,optarg);
      break;

    case 'o':
      strcpy(outfile,optarg);
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

    default:
      return 0;
    }
  }

  // Read data
  s=read_spectrogram(prefix,isub,nsub,f0,df0,nbin);
  printf("Read\n");
  // Write data
  write_spectrogram(s,outfile);
  printf("Written\n");

  // Free
  //  free(s.z);
  //  free(s.mjd);
  //  free(s.length);

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

// Compute Date from Julian Day
void mjd2nfd(double mjd,char *nfd)
{
  double f,jd,dday;
  int z,alpha,a,b,c,d,e;
  int year,month,day,hour,min;
  float sec,x;

  jd=mjd+2400000.5;
  jd+=0.5;

  z=floor(jd);
  f=fmod(jd,1.);

  if (z<2299161)
    a=z;
  else {
    alpha=floor((z-1867216.25)/36524.25);
    a=z+1+alpha-floor(alpha/4.);
  }
  b=a+1524;
  c=floor((b-122.1)/365.25);
  d=floor(365.25*c);
  e=floor((b-d)/30.6001);

  dday=b-d-floor(30.6001*e)+f;
  if (e<14)
    month=e-1;
  else
    month=e-13;

  if (month>2)
    year=c-4716;
  else
    year=c-4715;

  day=(int) floor(dday);
  x=24.0*(dday-day);
  x=3600.*fabs(x);
  sec=fmod(x,60.);
  x=(x-sec)/60.;
  min=fmod(x,60.);
  x=(x-min)/60.;
  hour=x;
  sec=floor(1000.0*sec)/1000.0;

  sprintf(nfd,"%04d-%02d-%02dT%02d:%02d:%06.3f",year,month,day,hour,min,sec);

  return;
}
