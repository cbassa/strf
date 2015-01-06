#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sgdp4h.h"
#include "satutl.h"
#include "rftime.h"

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
struct site get_site(int site_id);
void obspos_xyz(double mjd,double lng,double lat,float alt,xyz_t *pos,xyz_t *vel);
double gmst(double mjd);
double dgmst(double mjd);
double modulo(double x,double y);
void equatorial2horizontal(double mjd,struct site s,double ra,double de,double *azi,double *alt);

int main(int argc,char *argv[])
{
  int i,j,n=86400;
  int imode;
  double *mjd,mjd0=57028.0;
  struct point *p;
  struct site s;
  FILE *file;
  orbit_t orb;
  double dx,dy,dz,dvx,dvy,dvz,za;
  double r0,v0,ra0,de0,azi0,alt0;
  xyz_t satpos,satvel;
  double freq,freq0=2272.5e6;
  char nfd[32];
  double t=0.0,a,f;

  // Arrays
  mjd=(double *) malloc(sizeof(double)*n);
  p=(struct point *) malloc(sizeof(struct point)*n);
  
  // Get sites
  s=get_site(4171);

  // MJD range
  for (i=0;i<n;i++)
    mjd[i]=mjd0+(double) i/86400.0;

  // Get positions
  for (i=0;i<n;i++) 
    obspos_xyz(mjd[i],s.lng,s.lat,s.alt,&p[i].obspos,&p[i].obsvel);

  // Loop over objects
  file=fopen("/home/bassa/code/c/satellite/sattools/tle/catalog.tle","r");
  while (read_twoline(file,19822,&orb)==0) {
    // Initialize
    imode=init_sgdp4(&orb);
    if (imode==SGDP4_ERROR)
      printf("Error\n");
      
    // Loop over points
    for (i=0;i<n;i++) {
      // Get satellite position
      satpos_xyz(mjd[i]+2400000.5,&satpos,&satvel);

      // Observer position
      dx=satpos.x-p[i].obspos.x;  
      dy=satpos.y-p[i].obspos.y;
      dz=satpos.z-p[i].obspos.z;
      dvx=satvel.x-p[i].obsvel.x;
      dvy=satvel.y-p[i].obsvel.y;
      dvz=satvel.z-p[i].obsvel.z;
      r0=sqrt(dx*dx+dy*dy+dz*dz);
      v0=(dvx*dx+dvy*dy+dvz*dz)/r0;
      ra0=modulo(atan2(dy,dx)*R2D,360.0);
      de0=asin(dz/r0)*R2D;
      za=acos((p[i].obspos.x*dx+p[i].obspos.y*dy+p[i].obspos.z*dz)/(r0*XKMPER))*R2D;

      if (za<90.0)
	printf("%14.8lf %f %f %10.0lf %f\n",mjd[i],r0,v0,(1.0-v0/C)*freq0,za);
    }
  }
  fclose(file);

  // Free
  free(mjd);
  free(p);

  return 0;
}

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

// Convert equatorial into horizontal coordinates
void equatorial2horizontal(double mjd,struct site s,double ra,double de,double *azi,double *alt)
{
  double h;

  h=gmst(mjd)+s.lng-ra;
  
  *azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(s.lat*D2R)-tan(de*D2R)*cos(s.lat*D2R))*R2D,360.0);
  *alt=asin(sin(s.lat*D2R)*sin(de*D2R)+cos(s.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

  return;
}

