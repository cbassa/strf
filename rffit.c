#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <cpgplot.h>
#include <getopt.h>
#include "sgdp4h.h"

#define LIM 80
#define NMAX 1024
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.137 // km
#define KG 0.07436680 // earth radii, earth masses, minutes
#define C 299792.458 // km/s
#define FLAT (1.0/298.257)

struct site {
  int id;
  double lat,lng;
  float alt;
  char observer[64];
} site;
struct point {
  char timestamp[24];
  double mjd,freq,v,freq0;
  float t,f,res;
  float flux;
  int flag,site_id,rsite_id;
  struct site s,r;
};
struct data {
  int n;
  struct point *p;
  int fitfreq;
  double mjdmin,mjdmax,mjd0;
  double freqmin,freqmax,fluxmin,fluxmax,f0,ffit;
  char satname[LIM];
} d;
orbit_t orb;
int fgetline(FILE *file,char *s,int lim);
struct data read_data(char *filename,int graves,float offset);
double date2mjd(int year,int month,double day);
void mjd2nfd(double mjd,char *nfd);
struct point decode_line(char *line);
double modulo(double,double);
double gmst(double);
double dgmst(double);
void obspos_xyz(double,struct site,xyz_t *,xyz_t *);
int velocity(orbit_t orb,double mjd,struct site s,double *v,double *azi,double *alt);
double altitude(orbit_t orb,double mjd,struct site s);
void deselect_inside(float x0,float y0,float x,float y);
void highlight(float x0,float y0,float x,float y,int flag);
void deselect_outside(float xmin,float ymin,float xmax,float ymax);
void deselect_nearest(float x,float y,float xmin,float ymin,float xmax,float ymax);
void save_data(float xmin,float ymin,float xmax,float ymax,char *filename);
void equatorial2horizontal(double mjd,struct site s,double ra,double de,double *azi,double *alt);
double chisq(double a[]);
void versafit(int m,int n,double *a,double *da,double (*func)(double *),double dchisq,double tol,char *opt);
double compute_rms(void);
void mjd2date(double mjd,int *year,int *month,double *day);
void print_tle(orbit_t orb,char *filename);
void search(void);
double fit_curve(orbit_t orb,int *ia);
double mjd2doy(double mjd,int *yr);
double doy2mjd(int year,double doy);
int identify_satellite_from_visibility(char *catalog,double altmin);

// Get observing site
struct site get_site(int site_id)
{
  int i=0,status;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64];
  struct site s;
  char *env,filename[LIM];

  env=getenv("ST_DATADIR");
  if(env==NULL||strlen(env)==0)
    env=".";
  sprintf(filename,"%s/data/sites.txt",env);

  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not found!\n");
    exit(0);
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    status=sscanf(line,"%4d %2s %lf %lf %f",
	   &id,abbrev,&lat,&lng,&alt);
    strcpy(observer,line+38);

    // Change to km
    alt/=1000.0;
    
    // Copy site
    if (id==site_id) {
      s.lat=lat;
      s.lng=lng;
      s.alt=alt;
      s.id=id;
      strcpy(s.observer,observer);
    }

  }
  fclose(file);

  return s;
}

// Select diagonal
void diagonal_select(float x0,float y0,float x1,float y1,int flag)
{
  int i;
  float v;
  float ymin,ymax;
  printf("%f %f %f %f\n",x0,y0,x1,y1);

  for (i=0;i<d.n;i++) {
    v=(d.p[i].t-x0)/(x1-x0);
    ymin=y0+v*(y1-y0)-3.0;
    ymax=y0+v*(y1-y0)+3.0;
    if (v>0.0 && v<1.0 && d.p[i].f>ymin && d.p[i].f<ymax && d.p[i].flag!=0)
      d.p[i].flag=flag;
  }

  return;
}

void format_tle(orbit_t orb,char *line1,char *line2)
{
  int i,csum;
  char sbstar[]=" 00000-0",bstar[13];
  char csumstr[2];

  // Format Bstar term
  if (fabs(orb.bstar)>1e-9) {
    sprintf(bstar,"%11.4e",10*orb.bstar);
    sbstar[0] = bstar[0];  sbstar[1] = bstar[1];  sbstar[2] = bstar[3];  sbstar[3] = bstar[4];
    sbstar[4] = bstar[5];  sbstar[5] = bstar[6];  sbstar[6] = bstar[8];  sbstar[7] = bstar[10];  sbstar[8] = '\0';
  }
  // Print lines
  sprintf(line1,"1 %05dU          %2d%012.8f  .00000000  00000-0 %8s 0    0",orb.satno,orb.ep_year-2000,orb.ep_day,sbstar);
  sprintf(line2,"2 %05d %8.4f %8.4f %07.0f %8.4f %8.4f %11.8f    0",orb.satno,DEG(orb.eqinc),DEG(orb.ascn),1E7*orb.ecc,DEG(orb.argp),DEG(orb.mnan),orb.rev);

  // Compute checksums
  for (i=0,csum=0;i<strlen(line1);i++) {
    if (isdigit(line1[i]))
      csum+=line1[i]-'0';
    else if (line1[i]=='-')
      csum++;
  }
  sprintf(csumstr,"%d",csum%10);
  strcat(line1,csumstr);
  for (i=0,csum=0;i<strlen(line2);i++) {
    if (isdigit(line2[i]))
      csum+=line2[i]-'0';
    else if (line2[i]=='-')
      csum++;
  }
  sprintf(csumstr,"%d",csum%10);
  strcat(line2,csumstr);

  return;
}


int identify_satellite_from_doppler(char *catalog,double rmsmax)
{
  int i=0,flag=0;
  FILE *fp;
  double rms,rmsmin;
  int ia[]={0,0,0,0,0,0};
  int satno=0,imode;
  double v,alt,azi,mjdmid;

  // Open catalog
  fp=fopen(catalog,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",catalog);

  // Loop over TLEs
  while (read_twoline(fp,0,&orb)==0) {
    // Initialize
    imode=init_sgdp4(&orb);
    if (imode==SGDP4_ERROR)
      printf("Error\n");

    //    velocity(orb,d.p[d.n/2].mjd,d.p[d.n/2].s,&v,&azi,&alt);
    //    if (alt<0.0)
    //      printf("Continue?\n");
    //      continue;
    rms=fit_curve(orb,ia);
    if (flag==0 || rms<rmsmin) {
      rmsmin=rms;
      satno=orb.satno;
      flag=1;
    }
    if (rms<rmsmax) {
      printf("%05d %.3f kHz %.6f MHz\n",orb.satno,rms,d.ffit/1000.0);
      i++;
    }
  }
  rewind(fp);

  // Plot results
  if (i>0) {
    printf("Identified %d candidate(s), best fitting satellite is %05d.\n",i,satno);
    read_twoline(fp,satno,&orb);
    rms=fit_curve(orb,ia);
  } else {
    printf("No candidates found.\n");
    satno=-1;
  }
  fclose(fp);

  return satno;
}

int identify_satellite_from_visibility(char *catalog,double altmin)
{
  int i=0,flag=0,nalt,nsel;
  FILE *fp;
  int satno=0,imode;
  double alt,frac;

  // Open catalog
  fp=fopen(catalog,"rb");
  if (fp==NULL)
    fatal_error("File open failed for reading %s\n",catalog);

  // Loop over TLEs
  while (read_twoline(fp,0,&orb)==0) {
    // Initialize
    imode=init_sgdp4(&orb);
    if (imode==SGDP4_ERROR)
      printf("Error\n");

    if (orb.rev<10.0)
      continue;
    
    // Loop over observations
    for (i=0,nalt=0,nsel=0;i<d.n;i++) {
      if (d.p[i].flag==2) {
	alt=altitude(orb,d.p[i].mjd,d.p[i].s);
	nsel++;
	if (alt>altmin)
	  nalt++;
      }
    }
    frac=(float) nalt/(float) nsel;
    if (nsel==nalt)
      printf("%5d %d/%d %.4f\n",orb.satno,nalt,nsel,frac);
  }
  rewind(fp);

  fclose(fp);

  return satno;
}

void usage()
{
  printf("dpplot -d <data file> -c [tle catalog] -i [satno] -h\n\ndata file:    Tabulated doppler curve\ntle catalog:  Catalog with TLE's (optional)\nsatno:        Satellite to load from TLE catalog (optional)\n\n");

  printf("rffit: fit RF observations\n\n");
  printf("-d <datafile>   Input data file with RF measurements\n");
  printf("-c <catalog>    Catalog with TLE's [optional]\n");
  printf("-i <satno>      NORAD ID of satellite to load\n");
  printf("-s <site>       Site ID\n");
  printf("-g              GRAVES data\n");
  printf("-m <offset>     Frequency offset to apply [Hz]\n");
  printf("-h              This help\n");
  
  return;
}

int main(int argc,char *argv[])
{
  int i,j,flag,redraw=1,plot_curve=1,plot_type=1,residuals=0,elset=0;
  int imode,year,style,color;
  char xlabel[64],ylabel[32],text[64],freqlist[128];
  int site_id=4171;
  float xmin,xmax,ymin,ymax;
  float xminsel,xmaxsel,yminsel,ymaxsel;
  float x0,y0,x,y;
  double mjd,v,v1,azi,alt,rms=0.0,day,mjdtca=56658.0,altmin=0.0;
  float t,f,vtca,foffset=0.0;
  char c,nfd[32]="2014-01-01T00:00:00";
  int mode=0,posn=0,click=0;
  char *catalog,*datafile,filename[64],string[64],bstar[10]=" 00000-0";
  int arg=0,nobs=0;
  FILE *fp,*std,*fpres;
  char line0[72],line1[72],line2[72];
  int ia[]={0,0,0,0,0,0,0};
  float dx[]={0.1,0.1,0.35,0.35,0.6,0.6,0.85},dy[]={0.0,-0.25,0.0,-0.25,0.0,-0.25,0.0};
  int satno=-1,status;
  struct site s0,s1;
  int site_number[16],nsite=0,graves=0;
  char *env;

  // Get site
  env=getenv("ST_COSPAR");
  if (env!=NULL) {
    site_id=atoi(env);
  } else {
    printf("ST_COSPAR environment variable not found.\n");
  }
  
  env=getenv("ST_DATADIR");
  if(env==NULL||strlen(env)==0)
    env=".";
  // Decode options
  while ((arg=getopt(argc,argv,"d:c:i:hs:gm:"))!=-1) {
    switch(arg) {
    case 'd':
      datafile=optarg;
      break;

    case 'c':
      catalog=optarg;
      break;

    case 'i':
      satno=atoi(optarg);
      break;

    case 'h':
      usage();
      return 0;
      break;

    case 'm':
      foffset=atof(optarg)/1000.0;
      break;
      
    case 's':
      site_id=atoi(optarg);
      break;

    case 'g':
      graves=1;
      break;

    default:
      usage();
      return 0;
    }
  }
  
  // Read data
  d=read_data(datafile,graves,foffset);
  d.fitfreq=1;

  // Set graves frequency
  if (graves==1) {
    d.fitfreq=0;
    d.ffit=143050.0;
  }

  // Count number of sites and assign colors
  for (i=0;i<d.n;i++) {
    // Check if site is assigned
    for (j=0,flag=0;j<nsite;j++) 
      if (site_number[j]==d.p[i].site_id)
	flag=1;
    
    // Not assigned
    if (flag==0) {
      site_number[nsite]=d.p[i].site_id;
      nsite++;
    }
    if (nsite>=16) {
      printf("Too many observing sites.\n");
      return 0;
    }
  }

  // Set default observing site
  site=get_site(site_id);

  // Read TLE
  if (satno>=0) {
    fp=fopen(catalog,"rb");
    if (fp==NULL)
      fatal_error("File open failed for reading %s\n",catalog);
    status=read_twoline(fp,satno,&orb);
    if (status==-1) {
      printf("No elements found for %5d\n",satno);
      satno=-1;
    } 
    fclose(fp);
  }

  if (freopen("/tmp/stderr.txt","w",stderr)==NULL)
    fprintf(stderr,"Failed to redirect stderr\n");

  cpgopen("/xs");
  cpgask(0);

  xmin=d.mjdmin-d.mjd0-0.005;
  xmax=d.mjdmax-d.mjd0+0.005;
  if (graves==0) {
    ymin=-12.0*d.f0/C;
    ymax=12*d.f0/C;
  } else if (graves==1) {
    ymin=-20.0*d.f0/C;
    ymax=20*d.f0/C;
  }
    
  //  ymin=d.freqmin;
  //  ymax=d.freqmax;

  day=mjd2doy(d.mjd0,&year);
  sprintf(xlabel,"MJD - %5.0f (%02d%03.0f)",d.mjd0,year-2000,floor(day));
  sprintf(ylabel,"Frequency - %.0f kHz",d.f0);

  // For ever loop
  for (;;) {
    if (redraw==1) {
      // Plot buttons
      cpgpage();
      cpgsvp(0.1,0.95,0.0,0.2);
      cpgswin(0.0,1.0,-0.5,0.5);

      // Buttons
      cpgtext(0.12,-0.05,"Inclination");
      cpgtext(0.372,-0.05,"Eccentricity");
      cpgtext(0.62,-0.05,"Mean Anomaly");
      cpgtext(0.87,-0.05,"B\\u*\\d");
      cpgtext(0.12,-0.3,"Ascending Node");
      cpgtext(0.37,-0.3,"Arg. of Perigee");
      cpgtext(0.62,-0.3,"Mean Motion");

      // Toggles
      for (i=0;i<7;i++) {
	cpgpt1(dx[i],dy[i],19);
	if (ia[i]==1) {
	  cpgsci(2);
	  cpgpt1(dx[i],dy[i],16);
	  cpgsci(1);
	}
      }

      // Sky plot
      cpgsvp(0.715,0.95,0.5,0.99);
      cpgwnad(-90.0,90.0,-90.0,90.0);
      cpgbox("BC",0.,0,"BC",0.,0);
      cpgsfs(2);
      cpgcirc(0.0,0.0,90.0);
      cpgpt1(0.0,0.0,2);

      // Plot orbit
      if (satno>0 && plot_curve==1) {
	// Initialize
	imode=init_sgdp4(&orb);
	if (imode==SGDP4_ERROR)
	  break;
	
	cpgsci(15);
	//	for (mjd=d.mjd0,i=0;mjd<d.mjd0+xmax;mjd+=1.0/1440.0,i++) {
	for (i=0;i<NMAX;i++) {
	  mjd=xmin+d.mjd0+(xmax-xmin)*(float) i/(float) (NMAX-1);
	  velocity(orb,mjd,site,&v,&azi,&alt);

	  // Get TCA
	  if (i>0) {
	    if (vtca*v<0.0) {
	      mjdtca=mjd;
	    }
	  }

	  x=(90-alt)*sin(azi*D2R);
	  y=-(90-alt)*cos(azi*D2R);
	  if (i==0 || alt<0.0)
	    cpgmove(x,y);
	  else
	    cpgdraw(x,y);
	  vtca=v;
	}
	mjd2nfd(mjdtca,nfd);
	cpgsci(1);
	cpgsls(1);

	// Plot points
	for (i=0;i<d.n;i++) {
	  velocity(orb,d.p[i].mjd,site,&v,&azi,&alt);
	  x=(90-alt)*sin(azi*D2R);
	  y=-(90-alt)*cos(azi*D2R);
	  if (alt<0.0)
	    continue;
	  if (d.p[i].flag==1) {
	    cpgsci(1);
	    if (d.p[i].rsite_id==0)
	      cpgpt1(x,y,17);
	    else
	      cpgpt1(x,y,2);
	  } else if (d.p[i].flag==2) {
	    cpgsci(1);
	    if (d.p[i].rsite_id==0)
	      cpgpt1(x,y,17);
	    else
	      cpgpt1(x,y,2);
	    cpgsci(2);
	    cpgpt1(x,y,4);
	  }
	}
	cpgsci(1);
      }

      // Diagnostics
      cpgsvp(0.715,0.95,0.2,0.5);
      cpgwnad(0.0,1.0,0.0,1.0);
      
      sprintf(text,"Measurements: %d",nobs);
      cpgtext(0.0,1.0,text);
      sprintf(text,"Frequency: %.3f MHz",d.ffit/1000.0);
      cpgtext(0.0,0.85,text);
      sprintf(text,"rms: %.3f kHz",rms);
      cpgtext(0.0,0.7,text);
      sprintf(text,"TCA: %s",nfd);
      cpgtext(0.0,0.55,text);
      sprintf(text,"%s (%04d)",site.observer,site.id);
      cpgtext(0.0,0.4,text);

      // Plot site numbers
      for (j=0;j<nsite;j++) {
	sprintf(text,"%04d",site_number[j]);
	cpgsci(j+2);
	if (j<5)
	  cpgtext(0.25*j,0.25,text);
	else if (j<10)
	  cpgtext(0.25*(j-5),0.15,text);
	else if (j<15)
	  cpgtext(0.25*(j-10),0.05,text);
      }
      cpgsci(1);

      // Initialize
      cpgsvp(0.1,0.65,0.2,0.9);
      cpgswin(xmin,xmax,ymin/d.f0*C,ymax/d.f0*C);
      cpgbox("",0.,0,"CTSM",0.,0);
      cpgmtxt("R",2.5,0.5,0.5,"Velocity (km/s)");
      cpgswin(xmin,xmax,ymin,ymax);
      cpgbox("BCTSN",0.,0,"BTSN",0.,0);
      //      cpgenv(xmin,xmax,ymin,ymax,0,0);
      cpglab(xlabel,ylabel,"");

      // Plot orbit
      if (satno>0 && plot_curve==1 && residuals==0) {

	// Plot tle
	format_tle(orb,line1,line2);
	cpgmtxt("T",2.0,0.0,0.0,line1);
	cpgmtxt("T",1.0,0.0,0.0,line2);

	// Initialize
	imode=init_sgdp4(&orb);
	if (imode==SGDP4_ERROR)
	  break;
	
	// Loop over sites for plotting model
	for (j=0;j<nsite;j++) {
	  s0=get_site(site_number[j]);
	  s1=get_site(d.p[0].rsite_id);
	  color=j+2;

	  for (i=0;i<NMAX;i++) {
	    mjd=xmin+d.mjd0+(xmax-xmin)*(float) i/(float) (NMAX-1);
	    t=(float) (mjd-d.mjd0);
	    if (d.p[0].rsite_id!=0) {
	      velocity(orb,mjd,s1,&v1,&azi,&alt);
	      velocity(orb,mjd,s0,&v,&azi,&alt);
	      f=(float) ((1.0-v/C)*(1.0-v1/C)*d.ffit-d.f0);
	    } else {
	      velocity(orb,mjd,s0,&v,&azi,&alt);
	      f=(float) ((1.0-v/C)*d.ffit-d.f0);
	    }
	    
	    if (alt>0.0) {
	      cpgsls(1);
	      cpgsci(color);
	    } else {
	      cpgsls(2);
	      cpgsci(14);
	    }
	    
	    if (i==0) 
	      cpgmove(t,f);
	    else
	      cpgdraw(t,f);
	  }
	  cpgsci(1);
	  cpgsls(1);
	}
      } 
      // Plot selected points
      for (i=0;i<d.n;i++) {
	for (j=0;j<nsite;j++) 
	  if (d.p[i].site_id==site_number[j])
	    break;
	
	color=j+2;
	style=17;
	
	x=d.p[i].t;
	y=d.p[i].f;
	if (d.p[i].flag==1) {
	  cpgsci(color);
	  if (d.p[i].rsite_id==0)
	    cpgpt1(x,y,style);
	  else
	    cpgpt1(x,y,2);
	} else if (d.p[i].flag==2) {
	  cpgsci(color);
	  if (d.p[i].rsite_id==0)
	    cpgpt1(x,y,style);
	  else
	    cpgpt1(x,y,2);
	  cpgsci(1);
	  cpgpt1(x,y,4);
	}
      }
      cpgsci(1);
      redraw=0;
    }
    
    // Get cursor
    cpgband(mode,posn,x0,y0,&x,&y,&c);

    // Quit
    if (c=='q' || c=='Q')
      break;
   
    // Toggle curve
    if (c=='p') {
      plot_curve=(plot_curve==1) ? 0 : 1;
      redraw=1;
    }

    // Toggles
    if (isdigit(c) && c-'0'>=1 && c-'0'<8) {
      if (ia[c-49]==0) 
	ia[c-49]=1;
      else if (ia[c-49]==1) 
	ia[c-49]=0;
      redraw=1;
    }

    // Change
    if (c=='c') {
      printf("(1) Inclination,     (2) Ascending Node,   (3) Eccentricity,\n(4) Arg. of Perigee, (5) Mean Anomaly,     (6) Mean Motion,\n(7) B* drag,         (8) Epoch,            (9) Satellite ID\n\nWhich parameter to change: ");
      status=scanf("%i",&i);
      if (i>=0 && i<=9) {
	printf("\nNew value: ");
	if (fgets(string,64,stdin)==NULL)
	  fprintf(stderr,"Failed to read string\n");
	status=scanf("%s",string);
	//	if (i==0) strcpy(d.satname,string);
	if (i==1) orb.eqinc=RAD(atof(string));
	if (i==2) orb.ascn=RAD(atof(string));
	if (i==3) orb.ecc=atof(string);
	if (i==4) orb.argp=RAD(atof(string));
	if (i==5) orb.mnan=RAD(atof(string));
	if (i==6) orb.rev=atof(string);
	if (i==7) orb.bstar=atof(string);
	if (i==8) {
	  orb.ep_year=2000+(int) floor(atof(string)/1000.0);
	  orb.ep_day=atof(string)-1000*floor(atof(string)/1000.0);
	}
	if (i==9) orb.satno=atoi(string);
	if (i==0) {
	  printf("%f\n",d.ffit);
	  d.ffit=atof(string);
	  d.fitfreq=0;
	}
	redraw=1;
      }
      printf("\n================================================================================\n");
    }

    // Move
    if (c=='m') {
      printf("Provide frequency offset (kHz): ");
      status=scanf("%lf",&v);

      // Loop over points
      for (i=0;i<d.n;i++) {
	if (d.p[i].flag==2) {
	  d.p[i].f+=v;
	  d.p[i].freq+=v;
	}
      }
      click=0;
      redraw=1;
      printf("\n================================================================================\n");
      continue;
    }

    // Flux limit
    if (c=='l') {
      env=getenv("ST_DATADIR");
      if(env==NULL||strlen(env)==0)
        env=".";
      sprintf(freqlist,"%s/data/frequencies.txt",env);  
      fp=fopen(freqlist,"a");
      fprintf(fp,"%05d %lf\n",orb.satno,d.ffit/1000.0);
      fclose(fp);
      printf("Logged %05d at %lf to %s\n",orb.satno,d.ffit/1000.0,freqlist); 
      printf("\n================================================================================\n");
    }

    // Fit
    if (c=='f') {
      // Count points
      for (i=0,nobs=0;i<d.n;i++)
	if (d.p[i].flag==2)
	  nobs++;
      if (satno<0) {
	printf("No elements loaded!\n");
      } else if (nobs==0) {
	printf("No points selected!\n");
      } else {
	rms=fit_curve(orb,ia);
	//printf("rms: %.3f kHz, cf: %.6f MHz, TCA: %s\n",rms,d.ffit/1000.0,nfd);
	printf("%05d %.3f %.3f %s %04d %02d%010.6f\n",orb.satno,d.ffit/1000.0,rms,nfd,d.p[0].site_id,orb.ep_year-2000,orb.ep_day);
	redraw=1;
      }
    }

    // Get TLE
    if (c=='g') {
      printf("Get TLE from catalog, provide satellite number: ");
      status=scanf("%d",&satno);

      // Read TLE
      fp=fopen(catalog,"rb");
      if (fp==NULL)
	fatal_error("File open failed for reading %s\n",catalog);
      status=read_twoline(fp,satno,&orb);
      fclose(fp);
      if (status==-1) {
	printf("No elements found for %5d\n",satno);
	satno=-1;
      } else {
	print_orb(&orb);
	d.ffit=d.f0;
	redraw=1;
      }
      printf("\n================================================================================\n");
    }

    // Identify
    if (c=='i') {
      if (graves==0) {
	printf("rms limit (kHz): ");
	status=scanf("%lf",&rms);
      } else {
	printf("Using 0.1 kHz rms limit\n");
	rms=0.1;
      }
      satno=identify_satellite_from_doppler(catalog,rms);
      if (satno>0) {
	redraw=1;
	plot_curve=1;
      }
      printf("\n================================================================================\n");
    }

    // Identify
    if (c=='I') {
      printf("Above altitude (deg): ");
      status=scanf("%lf",&altmin);
      satno=identify_satellite_from_visibility(catalog,altmin);
      printf("\n================================================================================\n");
    }

    // Parameter search
    //    if (c=='S') {
    //      search();
    //      redraw=1;
    //    }

    // Write TLE
    if (c=='w') {
      printf("TLE filename to write: ");
      status=scanf("%s",filename);
      print_tle(orb,filename);
      printf("\n================================================================================\n");
    }

    // Reread tle
    if (c=='R') {
        // Read TLE
      fp=fopen(catalog,"rb");
      if (fp==NULL)
	fatal_error("File open failed for reading %s\n",catalog);
      read_twoline(fp,satno,&orb);
      print_orb(&orb);
      printf("\n================================================================================\n");
      fclose(fp);
      d.ffit=d.f0;
      redraw=1;
    }

    // Diagonal select
    if (c=='k') {
      printf("Selecting diagonal\n");
      diagonal_select(xmin,ymin,xmax,ymax,2);
      redraw=1;
    }

    // Toggle residuals
    if (c=='j') {
      fpres=fopen("residuals.dat","w");
      for (i=0;i<d.n;i++) 
	fprintf(fpres,"%14.8lf %lf %lf\n",d.p[i].mjd,d.p[i].freq,d.p[i].res);
      fclose(fpres);
      if (residuals==0)
	residuals=1;
      else if (residuals==1)
	residuals=0;
      redraw=1;
    }

    // Delete
    if (c=='X') {
      if (click==0) {
	deselect_nearest(x,y,xmin,ymin,xmax,ymax);

	click=0;
	redraw=1;
      } 
    }

    // Zoom
    if (c=='z') {
      click=1;
      mode=2;
    }
    // Delete box
    if (c=='d') {
      click=2;
      mode=2;
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
      } else if (click==2 && mode==2) {
	xminsel=(x0<x) ? x0 : x;
	xmaxsel=(x0>x) ? x0 : x;
	yminsel=(y0<y) ? y0 : y;
	ymaxsel=(y0>y) ? y0 : y;
	deselect_inside(xminsel,yminsel,xmaxsel,ymaxsel);
	click=0;
	mode=0;
	redraw=1;
      } else if (click==3 && mode==2) {
	xminsel=(x0<x) ? x0 : x;
	xmaxsel=(x0>x) ? x0 : x;
	yminsel=(y0<y) ? y0 : y;
	ymaxsel=(y0>y) ? y0 : y;
	printf("%f %f %f %f\n",xminsel,xmaxsel,yminsel,ymaxsel);
	click=0;
	mode=0;
	redraw=1;
      } else {
	click=0;
	mode=0;
	redraw=1;
      }
    }

    // Highlight
    if (c=='s') {
      highlight(xmin,ymin,xmax,ymax,2);
      for (i=0,nobs=0;i<d.n;i++)
	if (d.p[i].flag==2)
	  nobs++;
      click=0;
      mode=0;
      redraw=1;
    }

    // Highlight
    if (c=='H') {
      highlight(xmin,ymin,xmax,ymax,1);
      click=0;
      mode=0;
      redraw=1;
    }

    // Deselect all except highlighted
    if (c=='x') {
      deselect_inside(xmin,ymin,xmax,ymax);
      click=0;
      mode=0;
      redraw=1;
    }

    // Invert selection
    if (c=='T') {
      for (i=0;i<d.n;i++) {
	if (d.p[i].flag==2)
	  d.p[i].flag=1;
	else if (d.p[i].flag==1)
	  d.p[i].flag=2;
      }
      click=0;
      redraw=1;
    }

    // Deselect highlighted
    if (c=='D') {
      for (i=0;i<d.n;i++) {
	if (d.p[i].flag==2)
	  d.p[i].flag=0;
      }
      click=0;
      mode=0;
      redraw=1;
    }

    // Save
    if (c=='S') {
      printf("%s_%.3f_%04d_%05d.dat\n",nfd,d.ffit/1000.0,d.p[0].site_id,satno);
      printf("Save highlighted points, provide filename: ");
      status=scanf("%s",filename);
      save_data(xmin,ymin,xmax,ymax,filename);
      printf("\n================================================================================\n");
    }

    // Unselect
    if (c=='U') {
      for (i=0;i<d.n;i++)
	d.p[i].flag=1;
      redraw=1;
    }

    // Unselect
    if (c=='u') {
      for (i=0;i<d.n;i++) 
	if (d.p[i].flag==2)
	  d.p[i].flag=1;
      redraw=1;
    }

    // Default tle
    if (c=='t') {
      orb.satno=99999;
      strcpy(orb.desig,"13900A");
      orb.ep_day=mjd2doy(0.5*(d.mjdmin+d.mjdmax),&orb.ep_year);
      satno=orb.satno;
      if (elset==0) {
	orb.eqinc=0.5*M_PI;
	orb.ascn=modulo(gmst(0.5*(d.mjdmin+d.mjdmax))+site.lng,360.0)*D2R;
	orb.ecc=0.0001;
	orb.argp=0.0;
	orb.mnan=site.lat*D2R;
	orb.rev=14.0;
	orb.bstar=0.5e-4;
	printf("LEO orbit\n");
      } else if (elset==1) {
	orb.eqinc=20.0*D2R;
	orb.ascn=0.0;
	orb.ecc=0.7;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=2.25;
	orb.bstar=0.0;
	printf("GTO orbit\n");
      } else if (elset==2) {
	orb.eqinc=10.0*D2R;
	orb.ascn=0.0;
	orb.ecc=0.0;
	orb.argp=0.0;
	orb.mnan=0.0;
	orb.rev=1.0027;
	orb.bstar=0.0;
	printf("GSO orbit\n");
      } else if (elset==3) {
	orb.eqinc=63.434*D2R;
	orb.ascn=0.0;
	orb.ecc=0.71;
	orb.argp=270.0*D2R;
	orb.mnan=0.0;
	orb.rev=2.006;
	orb.bstar=0.0;
	printf("HEO orbit\n");
      }
      elset++;
      if (elset>3)
	elset=0;
      print_orb(&orb);
      printf("\n================================================================================\n");
      click=0;
      redraw=1;
      continue;
    }
    /*
    // Default tle
    if (c=='t') {
      orb.satno=99999;
      orb.eqinc=0.5*M_PI;
      orb.ascn=0.0;
      orb.ecc=0.0;
      orb.argp=0.0;
      orb.mnan=0.0;
      orb.rev=14.0;
      orb.bstar=0.0;
      orb.ep_day=mjd2doy(0.5*(d.mjdmin+d.mjdmax),&orb.ep_year);
      satno=99999;
      print_orb(&orb);
      printf("\n================================================================================\n");
      click=0;
      redraw=1;
      continue;
    }
    */
    // Unzoom
    if (c=='r') {
      xmin=d.mjdmin-d.mjd0;
      xmax=d.mjdmax-d.mjd0;
      if (graves==0) {
	ymin=-12.0*d.f0/C;
	ymax=12*d.f0/C;
      } else if (graves==1) {
	ymin=-20.0*d.f0/C;
	ymax=20*d.f0/C;
      }
      mode=0;
      click=0;
      redraw=1;
      continue;
    }

    

    // Help
    if (c=='h') {
      printf("Usage:\n");
      printf("================================================================================\n");
      printf("q   Quit\n");
      printf("p   Toggle curve plotting\n");
      printf("1   Toggle fitting parameter (Inclination)\n");
      printf("2   Toggle fitting parameter (RA of ascending node)\n");
      printf("3   Toggle fitting parameter (Eccentricity)\n");
      printf("4   Toggle fitting parameter (Argyment of perigee)\n");
      printf("5   Toggle fitting parameter (Mean anomaly)\n");
      printf("6   Toggle fitting parameter (Mean motion)\n");
      printf("c   Change parameter\n");
      printf("m   Move highlighted points in frequency\n");
      printf("l   Select points on flux limit\n");
      printf("f   Fit highlighted points\n");
      printf("g   Get TLE from catalog\n");
      printf("i   Identify satellite from catalog based on Doppler curve\n");
      printf("I   Identify satellite from catalog based on visibility\n");
      printf("w   Write present TLE\n");
      printf("R   Reread TLE from catalog\n");
      printf("X   Delete nearest point (right mouse button)\n");
      printf("z   Start box to zoom\n");
      printf("d   Start box to delete points\n");
      printf("A   Zoom/delete points (left mouse button)\n");
      printf("h   Highlight points in present window\n");
      printf("x   Deselect all except highlighted\n");
      printf("T   Invert selection\n");
      printf("D   Delete highlighted points\n");
      printf("s   Save highlighted points into file\n");
      printf("U   Deselect all points\n");
      printf("u   Deselect highlighted points\n");
      printf("t   Load template tle\n");
      printf("r   Reset zoom\n");
      printf("h   This help\n");
      printf("================================================================================\n\n");
    }


    // Save
    x0=x;
    y0=y;
  }
  cpgclos();


  // Free
  free(d.p);

  fclose(stderr);

  return 0;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;

  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  if (c == '\t')
    c=' ';
  if (c == '\n')
    s[i++] = c;
  s[i] = '\0';
  return i;
}

// Decode line
struct point decode_line(char *line) 
{
  int year,month,iday,hour,min,sec,fsec;
  double day;
  struct point p;
  int site_id,rsite_id,nread,status;

  // Get timestamp
  if (line[5]=='.') {
    site_id=-1;
    rsite_id=-1;
    nread=sscanf(line,"%lf %lf %f %d %d",&p.mjd,&p.freq,&p.flux,&site_id,&rsite_id);
  } else {
    strncpy(p.timestamp,line,19);
    p.timestamp[19]='\0';

    // Get MJD, freq, flux
    site_id=-1;
    status=sscanf(line,"%2d/%2d/%2d %2d:%2d:%2d.%1d %lf %f %d",&year,&month,&iday,&hour,&min,&sec,&fsec,&p.freq,&p.flux,&site_id);

    day=iday+hour/24.0+min/1440.0+(sec+0.1*fsec)/86400.0;
    p.mjd=date2mjd(2000+year,month,day);
  }

  // Decode site
  p.s=get_site(site_id);
  p.site_id=site_id;
  if (rsite_id!=-1) {
    p.r=get_site(rsite_id);
    p.rsite_id=rsite_id;
  } else {
    p.rsite_id=0;
  }

  // Change to kHz
  p.freq*=1E-3;

  // Set flag
  p.flag=1;

  return p;
}

// Read data
struct data read_data(char *filename,int graves,float offset)
{
  int i=0;
  char line[LIM];
  FILE *file;
  struct data d;
  double c;
  float sum,ssum,w,df;

  // Open file
  file=fopen(filename,"r");
  if (file==NULL) {
    fprintf(stderr,"Failed to open %s\n",filename);
    exit(1);
  }

  // Count lines
  while (fgetline(file,line,LIM)>0) 
    i++;
  d.n=i;

  // Allocate
  d.p=(struct point *) malloc(sizeof(struct point)*d.n);

  // Rewind file
  rewind(file);

  // Read data
  i=0;
  while (fgetline(file,line,LIM)>0) 
    d.p[i++]=decode_line(line);

  // Close file
  fclose(file);

  // Add frequency offset
  for (i=0;i<d.n;i++)
    d.p[i].freq+=offset;

  d.mjdmin=d.mjdmax=d.p[0].mjd;
  d.freqmin=d.freqmax=d.p[0].freq;
  d.fluxmin=d.fluxmax=d.p[0].flux;
  for (i=1;i<d.n;i++) {
    if (d.p[i].flag>0) {
      if (d.p[i].mjd<d.mjdmin) d.mjdmin=d.p[i].mjd;
      if (d.p[i].mjd>d.mjdmax) d.mjdmax=d.p[i].mjd;
      if (d.p[i].freq<d.freqmin) d.freqmin=d.p[i].freq;
      if (d.p[i].freq>d.freqmax) d.freqmax=d.p[i].freq;
      if (d.p[i].flux<d.fluxmin) d.fluxmin=d.p[i].flux;
      if (d.p[i].flux>d.fluxmax) d.fluxmax=d.p[i].flux;
    }
  }
  c=0.1*(d.mjdmax-d.mjdmin);
  d.mjdmin-=c;
  d.mjdmax+=c;
  c=0.1*(d.freqmax-d.freqmin);
  d.freqmin-=c;
  d.freqmax+=c;
  c=0.02*(d.fluxmax-d.fluxmin);
  d.fluxmin-=c;
  d.fluxmax+=c;

  // Center time and frequency
  d.mjd0=floor(0.5*(d.mjdmax+d.mjdmin));
  d.f0=floor(0.5*(d.freqmax+d.freqmin));
  if (graves==1)
    d.f0=143050.0;

  // Compute times
  for (i=0;i<d.n;i++) {
    d.p[i].t=(float) (d.p[i].mjd-d.mjd0);
    d.p[i].f=(float) (d.p[i].freq-d.f0);
  }
  d.ffit=d.f0;

  return d;
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

// SGDP4 line-of-sight velocity
int velocity(orbit_t orb,double mjd,struct site s,double *v,double *azi,double *alt)
{
  double dx,dy,dz,dvx,dvy,dvz,r;
  double ra,de;
  xyz_t satpos,obspos,satvel,obsvel;

    
  // Loop over data points
  obspos_xyz(mjd,s,&obspos,&obsvel);
  satpos_xyz(mjd+2400000.5,&satpos,&satvel);

  dx=satpos.x-obspos.x;  
  dy=satpos.y-obspos.y;
  dz=satpos.z-obspos.z;
  dvx=satvel.x-obsvel.x;
  dvy=satvel.y-obsvel.y;
  dvz=satvel.z-obsvel.z;
  r=sqrt(dx*dx+dy*dy+dz*dz);
  *v=(dvx*dx+dvy*dy+dvz*dz)/r;

  ra=modulo(atan2(dy,dx)*R2D,360.0);
  de=asin(dz/r)*R2D;

  equatorial2horizontal(mjd,s,ra,de,azi,alt);

  return 0;
}

// SGDP4 algitude
double altitude(orbit_t orb,double mjd,struct site s)
{
  double dx,dy,dz,dvx,dvy,dvz,r;
  double ra,de,azi,alt;
  xyz_t satpos,obspos,satvel,obsvel;

  // Loop over data points
  obspos_xyz(mjd,s,&obspos,&obsvel);
  satpos_xyz(mjd+2400000.5,&satpos,&satvel);

  dx=satpos.x-obspos.x;  
  dy=satpos.y-obspos.y;
  dz=satpos.z-obspos.z;
  r=sqrt(dx*dx+dy*dy+dz*dz);

  ra=modulo(atan2(dy,dx)*R2D,360.0);
  de=asin(dz/r)*R2D;

  equatorial2horizontal(mjd,s,ra,de,&azi,&alt);

  return alt;
}

// Observer position
void obspos_xyz(double mjd,struct site s,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,sl,dtheta;

  sl=sin(s.lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*sl*sl);
  gc=1.0/ff+s.alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+s.alt/XKMPER;

  theta=gmst(mjd)+s.lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(s.lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(s.lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(s.lat*D2R)*XKMPER;
  vel->x=-gc*cos(s.lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(s.lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;
  
  return;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Deselect inside box
void deselect_inside(float x0,float y0,float x,float y)
{
  int i;
  float xmin,xmax,ymin,ymax;

  xmin=(x0<x) ? x0 : x;
  xmax=(x0>x) ? x0 : x;
  ymin=(y0<y) ? y0 : y;
  ymax=(y0>y) ? y0 : y;
  for (i=0;i<d.n;i++) 
    if (d.p[i].t>xmin && d.p[i].t<xmax && d.p[i].f>ymin && d.p[i].f<ymax && d.p[i].flag!=2)
      d.p[i].flag=0;

  return;
}


// Highlight
void highlight(float x0,float y0,float x,float y,int flag)
{
  int i;
  float xmin,xmax,ymin,ymax;

  xmin=(x0<x) ? x0 : x;
  xmax=(x0>x) ? x0 : x;
  ymin=(y0<y) ? y0 : y;
  ymax=(y0>y) ? y0 : y;
  for (i=0;i<d.n;i++) 
    if (d.p[i].t>xmin && d.p[i].t<xmax && d.p[i].f>ymin && d.p[i].f<ymax && d.p[i].flag!=0)
      d.p[i].flag=flag;

  return;
}

// Deselect outside box
void deselect_outside(float xmin,float ymin,float xmax,float ymax)
{
  int i;

  for (i=0;i<d.n;i++) 
    if (d.p[i].t<xmin || d.p[i].t>xmax || d.p[i].f<ymin || d.p[i].f>ymax)
      d.p[i].flag=0;

  return;
}

// Deselect nearest point
void deselect_nearest(float x,float y,float xmin,float ymin,float xmax,float ymax)
{
  int i,imin,flag;
  float r,rmin;
  float dx,dy;

  for (i=0,flag=0;i<d.n;i++) {
    if (d.p[i].flag>0) {
      dx=(x-d.p[i].t)/(xmax-xmin);
      dy=(y-d.p[i].f)/(ymax-ymin);;
      r=sqrt(dx*dx+dy*dy);
      if (flag==0) {
	imin=i;
	rmin=r;
	flag=1;
      } 
      if (r<rmin) {
	imin=i;
	rmin=r;
      }
    }
  }

  if (imin!=-1)
    d.p[imin].flag=0;

  return;
}

// Save data
void save_data(float xmin,float ymin,float xmax,float ymax,char *filename)
{
  int i,j;
  FILE *file;

  file=fopen(filename,"w");
  for (i=0,j=0;i<d.n;i++) {
    if (d.p[i].t>xmin && d.p[i].t<xmax && d.p[i].f>ymin && d.p[i].f<ymax && d.p[i].flag==2) {
      //      fprintf(file,"%s\t%14.3lf\t%8.3f\t%04d\n",d.p[i].timestamp,1000.0*d.p[i].freq,d.p[i].flux,d.p[i].site_id);
      if (d.p[i].rsite_id==0)
	fprintf(file,"%lf\t%14.3lf\t%8.3f\t%04d\n",d.p[i].mjd,1000.0*d.p[i].freq,d.p[i].flux,d.p[i].site_id);
      else
	fprintf(file,"%lf\t%14.3lf\t%8.3f\t%04d\t%04d\n",d.p[i].mjd,1000.0*d.p[i].freq,d.p[i].flux,d.p[i].site_id,d.p[i].rsite_id);
      j++;
    }
  }
  printf("%d points saved in %s\n",j,filename);
  fclose(file);

  return;
}

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,struct site s,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+s.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(s.lat*D2R)-tan(de*D2R)*cos(s.lat*D2R))*R2D,360.0);
  *alt=asin(sin(s.lat*D2R)*sin(de*D2R)+cos(s.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

  return;
}

// Chisq
double chisq(double a[])
{
  int i,imode;
  double *v,azi,alt,f,*v1,fac;
  double chisq;
  double sum1,sum2;
  
  // Allocate
  v=(double *) malloc(sizeof(double)*d.n);
  v1=(double *) malloc(sizeof(double)*d.n);

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day
  // a[6]: Bstar drag

  if (a[2]<0.0)
    a[2]=0.0;
  if (a[2]>=1.0)
    a[2]=0.999;
  if (a[5]<0.05)
    a[5]=0.05;

  // Set parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  // Initialize
  imode=init_sgdp4(&orb);
  if (imode==SGDP4_ERROR)
    printf("Error\n");

  // Loop over highlighted points
  for (i=0,sum1=0.0,sum2=0.0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      if (d.p[i].rsite_id!=0) {
	velocity(orb,d.p[i].mjd,d.p[i].r,&v1[i],&azi,&alt);
	velocity(orb,d.p[i].mjd,d.p[i].s,&v[i],&azi,&alt);
	fac=(1.0-v[i]/C)*(1.0-v1[i]/C);
	sum1+=fac*d.p[i].freq;
	sum2+=fac*fac;
      } else {
	velocity(orb,d.p[i].mjd,d.p[i].s,&v[i],&azi,&alt);
	fac=1.0-v[i]/C;
	sum1+=fac*d.p[i].freq;
	sum2+=fac*fac;
      }
    }
  }
  if (d.fitfreq==1)
    d.ffit=sum1/sum2;

  // Compute chisq
  for (i=0,chisq=0.0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      if (d.p[i].rsite_id!=0) 
	f=(1.0-v[i]/C)*(1.0-v1[i]/C)*d.ffit;
      else
	f=(1.0-v[i]/C)*d.ffit;
      chisq+=pow(d.p[i].freq-f,2);
    }
  }

  free(v);
  free(v1);

  return chisq;
}

// rms
double compute_rms(void)
{
  int i,imode,n;
  double v,v1,azi,alt,f;
  double rms;

  // Initialize
  imode=init_sgdp4(&orb);
  if (imode==SGDP4_ERROR)
    printf("Error\n");

  // Compute rms
  for (i=0,n=0,rms=0.0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      velocity(orb,d.p[i].mjd,d.p[i].s,&v,&azi,&alt);
      if (d.p[i].rsite_id!=0) {
	velocity(orb,d.p[i].mjd,d.p[i].r,&v1,&azi,&alt);
	f=(1.0-v/C)*(1.0-v1/C)*d.ffit;
      } else {
	f=(1.0-v/C)*d.ffit;
      }
      d.p[i].freq0=f;
      d.p[i].res=d.p[i].freq-f;
      rms+=pow(d.p[i].freq-f,2);

      n++;
    }
  }
  if (n>0)
    rms=sqrt(rms/(float) n);

  return rms;
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

  sprintf(nfd,"%04d-%02d-%02dT%02d:%02d:%02.0f",year,month,day,hour,min,sec);

  return;
}

// Compute Date from Julian Day
void mjd2date(double mjd,int *year,int *month,double *day)
{
  double f,jd;
  int z,alpha,a,b,c,d,e;

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

  *day=b-d-floor(30.6001*e)+f;
  if (e<14)
    *month=e-1;
  else
    *month=e-13;

  if (*month>2)
    *year=c-4716;
  else
    *year=c-4715;

  return;
}

// Print TLE
void print_tle(orbit_t orb,char *filename)
{
  int i,n;
  FILE *file;
  double mjdmin,mjdmax;
  int year,month;
  double day;
  char line1[70],line2[70];

  // Count number of points
  for (i=0,n=0;i<d.n;i++) {
    if (d.p[i].flag==2) {
      if (n==0) {
	mjdmin=d.p[i].mjd;
	mjdmax=d.p[i].mjd;
      }
      if (d.p[i].mjd<mjdmin) mjdmin=d.p[i].mjd;
      if (d.p[i].mjd>mjdmax) mjdmax=d.p[i].mjd;
      n++;
    }
  }

  // Write TLE
  file=fopen(filename,"w");
  format_tle(orb,line1,line2);
  fprintf(file,"%s\n%s\n",line1,line2);

  mjd2date(mjdmin,&year,&month,&day);
  fprintf(file,"# %4d%02d%05.2lf-",year,month,day);
  mjd2date(mjdmax,&year,&month,&day);
  fprintf(file,"%4d%02d%05.2lf, %d measurements, %.3lf kHz rms\n",year,month,day,n,compute_rms());
  fclose(file);

  return;
}

// Parameter search
void search(void) 
{
  int i,j,n;
  double a[7],da[7];
  FILE *file;
  double xmin,xmax,ymin,ymax;
  int nx,ny,status;
  
  // Get input
  printf("Provide mean motion estimate: ");
  status=scanf("%lf",&a[5]);
  printf("Provide inclination estimate: ");
  status=scanf("%lf",&a[0]);
  printf("Provide mean anomaly range, steps [min max nstep]: ");
  status=scanf("%lf %lf %d",&xmin,&xmax,&nx);
  printf("Provide ascending node range, steps [min max nstep]: ");
  status=scanf("%lf %lf %d",&ymin,&ymax,&ny);

  // Count highlighted points
  for (i=0,n=0;i<d.n;i++)
    if (d.p[i].flag==2)
      n++;
  
  // Loop
  printf("Starting parameter search\n");
  file=fopen("search.dat","w");
  for (i=0;i<nx;i++) {
    a[4]=xmin+(xmax-xmin)*(double) i/(double) (nx-1);
    for (j=0;j<ny;j++) {
      a[1]=ymin+(ymax-ymin)*(double) j/(double) (ny-1);

      //      a[1]=orb.ascn*R2D;
      a[2]=0.0;
      a[3]=0.0;
      //      a[4]=orb.mnan*R2D;
      da[0]=0.0;
      da[1]=0.0;
      da[2]=0.0;
      da[3]=0.0;
      da[4]=0.0;
      da[5]=0.0;
      da[6]=0.0;

      // Fit
      //versafit(n,6,a,da,chisq,0.0,1e-3,"n");

      orb.eqinc=RAD(a[0]);
      orb.ascn=RAD(modulo(a[1],360.0));
      orb.ecc=a[2];
      orb.argp=RAD(modulo(a[3],360.0));
      orb.mnan=RAD(modulo(a[4],360.0));
      orb.rev=a[5];
      orb.bstar=a[6];

      fprintf(file,"%8.4f %8.4f %f\n",a[4],a[1],compute_rms());
    }
    fprintf(file,"\n");
  }
  fclose(file);
  printf("Finished\n");

  return;
}

// Fit doppler curve
double fit_curve(orbit_t orb,int *ia)
{
  int i,n;
  double a[7],da[7];
  double db[7]={5.0,5.0,0.1,5.0,5.0,0.1,0.00001};
  double rms;

  a[0]=orb.eqinc*R2D;
  a[1]=orb.ascn*R2D;
  a[2]=orb.ecc;
  a[3]=orb.argp*R2D;
  a[4]=orb.mnan*R2D;
  a[5]=orb.rev;
  a[6]=orb.bstar;

  for (i=0;i<7;i++) {
    if (ia[i]==1)
      da[i]=db[i];
    else
      da[i]=0.0;
  }

  // Construct struct
  // a[0]: inclination
  // a[1]: RA of ascending node
  // a[2]: eccentricity
  // a[3]: argument of periastron
  // a[4]: mean anomaly
  // a[5]: revs per day
  // a[6]: Bstar drag

  // Count highlighted points
  for (i=0,n=0;i<d.n;i++)
    if (d.p[i].flag==2)
      n++;

  if (n>0)
    versafit(n,7,a,da,chisq,0.0,1e-5,"n");

  // Return parameters
  orb.eqinc=RAD(a[0]);
  orb.ascn=RAD(modulo(a[1],360.0));
  orb.ecc=a[2];
  orb.argp=RAD(modulo(a[3],360.0));
  orb.mnan=RAD(modulo(a[4],360.0));
  orb.rev=a[5];
  orb.bstar=a[6];

  //  print_tle(orb);
  rms=compute_rms();

  return rms;
}

// MJD to DOY
double mjd2doy(double mjd,int *yr)
{
  int year,month,k=2;
  double day,doy;
  
  mjd2date(mjd,&year,&month,&day);

  if (year%4==0 && year%400!=0)
    k=1;

  doy=floor(275.0*month/9.0)-k*floor((month+9.0)/12.0)+day-30;

  *yr=year;

  return doy;
}

// DOY to MJD
double doy2mjd(int year,double doy)
{
  int month,k=2;
  double day;

  if (year%4==0 && year%400!=0)
    k=1;

  month=floor(9.0*(k+doy)/275.0+0.98);
  
  if (doy<32)
    month=1;

  day=doy-floor(275.0*month/9.0)+k*floor((month+9.0)/12.0)+30.0;

  return date2mjd(year,month,day);
}
