#include <stdio.h>
#include <math.h>
#include "rftime.h"

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
