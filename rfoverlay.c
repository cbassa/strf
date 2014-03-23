#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sgdp4h.h"
#include "satutl.h"
#include "cpgplot.h"

#define LIM 80
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)
#define C 299792.458 // Speed of light in km/s

struct point {
  xyz_t obspos,obsvel;
};
struct site {
  int id;
  double lng,lat;
  float alt;
  char observer[64];
};

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Read a line of maximum length int lim from file FILE into string s
int fgetline(FILE *file,char *s,int lim)
{
  int c,i=0;
 
  while (--lim > 0 && (c=fgetc(file)) != EOF && c != '\n')
    s[i++] = c;
  //  if (c == '\n')
  //    s[i++] = c;
  s[i] = '\0';
  return i;
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

// Observer position
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,s,dtheta;

  s=sin(lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*s*s);
  gc=1.0/ff+alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+alt/XKMPER;

  theta=gmst(mjd)+lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(lat*D2R)*XKMPER;
  vel->x=-gc*cos(lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;

  return;
}

// Get observing site
struct site get_site(int site_id)
{
  int i=0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64];
  struct site s;
  char *env,filename[LIM];

  env=getenv("ST_DATADIR");
  sprintf(filename,"%s/data/sites.txt",env);

  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not found!\n");
    return;
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    sscanf(line,"%4d %2s %lf %lf %f",
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

// Plot overlay
void overlay(double *mjd,int n,int site_id)
{
  int i,imode,flag,satno,tflag;
  struct point *p;
  struct site s;
  FILE *file,*infile;
  orbit_t orb;
  xyz_t satpos,satvel;
  double dx,dy,dz,dvx,dvy,dvz,r,v,za;
  double freq,freq0;
  char line[LIM],text[8];

  // Get site
  s=get_site(site_id);

  // Allocate
  p=(struct point *) malloc(sizeof(struct point)*n);

  // Get observer position
  for (i=0;i<n;i++) 
    obspos_xyz(mjd[i],s.lng,s.lat,s.alt,&p[i].obspos,&p[i].obsvel);

  infile=fopen("frequencies.txt","r");
  while (fgetline(infile,line,LIM)>0) {
    sscanf(line,"%d %lf",&satno,&freq0);
    sprintf(text," %d",satno);
    // Loop over TLEs
    file=fopen("/home/bassa/code/c/satellite/sattools/tle/bulk.tle","r");
    while (read_twoline(file,satno,&orb)==0) {
      // Initialize
      imode=init_sgdp4(&orb);
      if (imode==SGDP4_ERROR)
	printf("Error\n");
      
      // Loop over points
      for (i=0,flag=0,tflag=0;i<n;i++) {
	// Get satellite position
	satpos_xyz(mjd[i]+2400000.5,&satpos,&satvel);
	
	dx=satpos.x-p[i].obspos.x;  
	dy=satpos.y-p[i].obspos.y;
	dz=satpos.z-p[i].obspos.z;
	dvx=satvel.x-p[i].obsvel.x;
	dvy=satvel.y-p[i].obsvel.y;
	dvz=satvel.z-p[i].obsvel.z;
	r=sqrt(dx*dx+dy*dy+dz*dz);
	v=(dvx*dx+dvy*dy+dvz*dz)/r;
	za=acos((p[i].obspos.x*dx+p[i].obspos.y*dy+p[i].obspos.z*dz)/(r*XKMPER))*R2D;
	
	freq=(1.0-v/C)*freq0;
	
	if (flag==0) {
	  cpgmove((float) i,(float) freq);
	  flag=1;
	} else {
	  cpgdraw((float) i,(float) freq);
	}
	
	if (za<90.0 && flag==0) {
	  flag=1;
	} else if (za>90.0 && flag==1) {
	  if (tflag==0) {
	    tflag=1;
	    cpgtext((float) i,(float) freq,text);
	  }
	  flag=0;
	}
      }
    }
    fclose(file);
  }
  fclose(infile);

  // Free
  free(p);

  return;
}
