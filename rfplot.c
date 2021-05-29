#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
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
  float sigma;
};
struct data {
  int n;
  double *x, *y, *ym;
} gd;

void dec2sex(double x,char *s,int f,int len);
void time_axis(double *mjd,int n,float xmin,float xmax,float ymin,float ymax);
void usage(void);
void plot_traces(struct trace *t,int nsat,float fcen);
struct trace fit_trace(struct spectrogram s,struct select sel,int site_id,int graves);
struct trace fit_gaussian_trace(struct spectrogram s,struct select sel,int site_id,int graves);
void convolve(float *y,int n,float *w,int m,float *z);
float gauss(float x,float w);
void quadfit(float x[],float y[],int n,float a[]);
struct trace locate_trace(struct spectrogram s,struct trace t,int site_id,float sigmamin,float wmin,float wmax,int graves);
void filter(struct spectrogram s,int site_id,float sigma,int graves);
void peakfind(struct spectrogram s,int site_id,int i0,int i1,int j0,int j1);
void versafit(int m,int n,double *a,double *da,double (*func)(double *),double dchisq,double tol,char *opt);
double chisq_gaussian(double a[]);
float fit_gaussian_point(struct spectrogram s,float x,float y,struct select sel,int site_id,int graves);

int main(int argc,char *argv[])
{
  struct spectrogram s;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};
  float cool_l[]={-0.5,0.0,0.17,0.33,0.50,0.67,0.83,1.0,1.7};
  float cool_r[]={0.0,0.0,0.0,0.0,0.6,1.0,1.0,1.0,1.0};
  float cool_g[]={0.0,0.0,0.0,1.0,1.0,1.0,0.6,0.0,1.0};
  float cool_b[]={0.0,0.3,0.8,1.0,0.3,0.0,0.0,0.0,1.0};
  float viridis_l[]={0.000000,0.003922,0.007843,0.011765,0.015686,0.019608,0.023529,0.027451,0.031373,0.035294,0.039216,0.043137,0.047059,0.050980,0.054902,0.058824,0.062745,0.066667,0.070588,0.074510,0.078431,0.082353,0.086275,0.090196,0.094118,0.098039,0.101961,0.105882,0.109804,0.113725,0.117647,0.121569,0.125490,0.129412,0.133333,0.137255,0.141176,0.145098,0.149020,0.152941,0.156863,0.160784,0.164706,0.168627,0.172549,0.176471,0.180392,0.184314,0.188235,0.192157,0.196078,0.200000,0.203922,0.207843,0.211765,0.215686,0.219608,0.223529,0.227451,0.231373,0.235294,0.239216,0.243137,0.247059,0.250980,0.254902,0.258824,0.262745,0.266667,0.270588,0.274510,0.278431,0.282353,0.286275,0.290196,0.294118,0.298039,0.301961,0.305882,0.309804,0.313725,0.317647,0.321569,0.325490,0.329412,0.333333,0.337255,0.341176,0.345098,0.349020,0.352941,0.356863,0.360784,0.364706,0.368627,0.372549,0.376471,0.380392,0.384314,0.388235,0.392157,0.396078,0.400000,0.403922,0.407843,0.411765,0.415686,0.419608,0.423529,0.427451,0.431373,0.435294,0.439216,0.443137,0.447059,0.450980,0.454902,0.458824,0.462745,0.466667,0.470588,0.474510,0.478431,0.482353,0.486275,0.490196,0.494118,0.498039,0.501961,0.505882,0.509804,0.513725,0.517647,0.521569,0.525490,0.529412,0.533333,0.537255,0.541176,0.545098,0.549020,0.552941,0.556863,0.560784,0.564706,0.568627,0.572549,0.576471,0.580392,0.584314,0.588235,0.592157,0.596078,0.600000,0.603922,0.607843,0.611765,0.615686,0.619608,0.623529,0.627451,0.631373,0.635294,0.639216,0.643137,0.647059,0.650980,0.654902,0.658824,0.662745,0.666667,0.670588,0.674510,0.678431,0.682353,0.686275,0.690196,0.694118,0.698039,0.701961,0.705882,0.709804,0.713725,0.717647,0.721569,0.725490,0.729412,0.733333,0.737255,0.741176,0.745098,0.749020,0.752941,0.756863,0.760784,0.764706,0.768627,0.772549,0.776471,0.780392,0.784314,0.788235,0.792157,0.796078,0.800000,0.803922,0.807843,0.811765,0.815686,0.819608,0.823529,0.827451,0.831373,0.835294,0.839216,0.843137,0.847059,0.850980,0.854902,0.858824,0.862745,0.866667,0.870588,0.874510,0.878431,0.882353,0.886275,0.890196,0.894118,0.898039,0.901961,0.905882,0.909804,0.913725,0.917647,0.921569,0.925490,0.929412,0.933333,0.937255,0.941176,0.945098,0.949020,0.952941,0.956863,0.960784,0.964706,0.968627,0.972549,0.976471,0.980392,0.984314,0.988235,0.992157,0.996078,1.000000};
  float viridis_r[]={0.267004,0.268510,0.269944,0.271305,0.272594,0.273809,0.274952,0.276022,0.277018,0.277941,0.278791,0.279566,0.280267,0.280894,0.281446,0.281924,0.282327,0.282656,0.282910,0.283091,0.283197,0.283229,0.283187,0.283072,0.282884,0.282623,0.282290,0.281887,0.281412,0.280868,0.280255,0.279574,0.278826,0.278012,0.277134,0.276194,0.275191,0.274128,0.273006,0.271828,0.270595,0.269308,0.267968,0.266580,0.265145,0.263663,0.262138,0.260571,0.258965,0.257322,0.255645,0.253935,0.252194,0.250425,0.248629,0.246811,0.244972,0.243113,0.241237,0.239346,0.237441,0.235526,0.233603,0.231674,0.229739,0.227802,0.225863,0.223925,0.221989,0.220057,0.218130,0.216210,0.214298,0.212395,0.210503,0.208623,0.206756,0.204903,0.203063,0.201239,0.199430,0.197636,0.195860,0.194100,0.192357,0.190631,0.188923,0.187231,0.185556,0.183898,0.182256,0.180629,0.179019,0.177423,0.175841,0.174274,0.172719,0.171176,0.169646,0.168126,0.166617,0.165117,0.163625,0.162142,0.160665,0.159194,0.157729,0.156270,0.154815,0.153364,0.151918,0.150476,0.149039,0.147607,0.146180,0.144759,0.143343,0.141935,0.140536,0.139147,0.137770,0.136408,0.135066,0.133743,0.132444,0.131172,0.129933,0.128729,0.127568,0.126453,0.125394,0.124395,0.123463,0.122606,0.121831,0.121148,0.120565,0.120092,0.119738,0.119512,0.119423,0.119483,0.119699,0.120081,0.120638,0.121380,0.122312,0.123444,0.124780,0.126326,0.128087,0.130067,0.132268,0.134692,0.137339,0.140210,0.143303,0.146616,0.150148,0.153894,0.157851,0.162016,0.166383,0.170948,0.175707,0.180653,0.185783,0.191090,0.196571,0.202219,0.208030,0.214000,0.220124,0.226397,0.232815,0.239374,0.246070,0.252899,0.259857,0.266941,0.274149,0.281477,0.288921,0.296479,0.304148,0.311925,0.319809,0.327796,0.335885,0.344074,0.352360,0.360741,0.369214,0.377779,0.386433,0.395174,0.404001,0.412913,0.421908,0.430983,0.440137,0.449368,0.458674,0.468053,0.477504,0.487026,0.496615,0.506271,0.515992,0.525776,0.535621,0.545524,0.555484,0.565498,0.575563,0.585678,0.595839,0.606045,0.616293,0.626579,0.636902,0.647257,0.657642,0.668054,0.678489,0.688944,0.699415,0.709898,0.720391,0.730889,0.741388,0.751884,0.762373,0.772852,0.783315,0.793760,0.804182,0.814576,0.824940,0.835270,0.845561,0.855810,0.866013,0.876168,0.886271,0.896320,0.906311,0.916242,0.926106,0.935904,0.945636,0.955300,0.964894,0.974417,0.983868,0.993248};
  float viridis_g[]={0.004874,0.009605,0.014625,0.019942,0.025563,0.031497,0.037752,0.044167,0.050344,0.056324,0.062145,0.067836,0.073417,0.078907,0.084320,0.089666,0.094955,0.100196,0.105393,0.110553,0.115680,0.120777,0.125848,0.130895,0.135920,0.140926,0.145912,0.150881,0.155834,0.160771,0.165693,0.170599,0.175490,0.180367,0.185228,0.190074,0.194905,0.199721,0.204520,0.209303,0.214069,0.218818,0.223549,0.228262,0.232956,0.237631,0.242286,0.246922,0.251537,0.256130,0.260703,0.265254,0.269783,0.274290,0.278775,0.283237,0.287675,0.292092,0.296485,0.300855,0.305202,0.309527,0.313828,0.318106,0.322361,0.326594,0.330805,0.334994,0.339161,0.343307,0.347432,0.351535,0.355619,0.359683,0.363727,0.367752,0.371758,0.375746,0.379716,0.383670,0.387607,0.391528,0.395433,0.399323,0.403199,0.407061,0.410910,0.414746,0.418570,0.422383,0.426184,0.429975,0.433756,0.437527,0.441290,0.445044,0.448791,0.452530,0.456262,0.459988,0.463708,0.467423,0.471133,0.474838,0.478540,0.482237,0.485932,0.489624,0.493313,0.497000,0.500685,0.504369,0.508051,0.511733,0.515413,0.519093,0.522773,0.526453,0.530132,0.533812,0.537492,0.541173,0.544853,0.548535,0.552216,0.555899,0.559582,0.563265,0.566949,0.570633,0.574318,0.578002,0.581687,0.585371,0.589055,0.592739,0.596422,0.600104,0.603785,0.607464,0.611141,0.614817,0.618490,0.622161,0.625828,0.629492,0.633153,0.636809,0.640461,0.644107,0.647749,0.651384,0.655014,0.658636,0.662252,0.665859,0.669459,0.673050,0.676631,0.680203,0.683765,0.687316,0.690856,0.694384,0.697900,0.701402,0.704891,0.708366,0.711827,0.715272,0.718701,0.722114,0.725509,0.728888,0.732247,0.735588,0.738910,0.742211,0.745492,0.748751,0.751988,0.755203,0.758394,0.761561,0.764704,0.767822,0.770914,0.773980,0.777018,0.780029,0.783011,0.785964,0.788888,0.791781,0.794644,0.797475,0.800275,0.803041,0.805774,0.808473,0.811138,0.813768,0.816363,0.818921,0.821444,0.823929,0.826376,0.828786,0.831158,0.833491,0.835785,0.838039,0.840254,0.842430,0.844566,0.846661,0.848717,0.850733,0.852709,0.854645,0.856542,0.858400,0.860219,0.861999,0.863742,0.865448,0.867117,0.868751,0.870350,0.871916,0.873449,0.874951,0.876424,0.877868,0.879285,0.880678,0.882046,0.883393,0.884720,0.886029,0.887322,0.888601,0.889868,0.891125,0.892374,0.893616,0.894855,0.896091,0.897330,0.898570,0.899815,0.901065,0.902323,0.903590,0.904867,0.906157};
  float viridis_b[]={0.329415,0.335427,0.341379,0.347269,0.353093,0.358853,0.364543,0.370164,0.375715,0.381191,0.386592,0.391917,0.397163,0.402329,0.407414,0.412415,0.417331,0.422160,0.426902,0.431554,0.436115,0.440584,0.444960,0.449241,0.453427,0.457517,0.461510,0.465405,0.469201,0.472899,0.476498,0.479997,0.483397,0.486697,0.489898,0.493001,0.496005,0.498911,0.501721,0.504434,0.507052,0.509577,0.512008,0.514349,0.516599,0.518762,0.520837,0.522828,0.524736,0.526563,0.528312,0.529983,0.531579,0.533103,0.534556,0.535941,0.537260,0.538516,0.539709,0.540844,0.541921,0.542944,0.543914,0.544834,0.545706,0.546532,0.547314,0.548053,0.548752,0.549413,0.550038,0.550627,0.551184,0.551710,0.552206,0.552675,0.553117,0.553533,0.553925,0.554294,0.554642,0.554969,0.555276,0.555565,0.555836,0.556089,0.556326,0.556547,0.556753,0.556944,0.557120,0.557282,0.557430,0.557565,0.557685,0.557792,0.557885,0.557965,0.558030,0.558082,0.558119,0.558141,0.558148,0.558140,0.558115,0.558073,0.558013,0.557936,0.557840,0.557724,0.557587,0.557430,0.557250,0.557049,0.556823,0.556572,0.556295,0.555991,0.555659,0.555298,0.554906,0.554483,0.554029,0.553541,0.553018,0.552459,0.551864,0.551229,0.550556,0.549841,0.549086,0.548287,0.547445,0.546557,0.545623,0.544641,0.543611,0.542530,0.541400,0.540218,0.538982,0.537692,0.536347,0.534946,0.533488,0.531973,0.530398,0.528763,0.527068,0.525311,0.523491,0.521608,0.519661,0.517649,0.515571,0.513427,0.511215,0.508936,0.506589,0.504172,0.501686,0.499129,0.496502,0.493803,0.491033,0.488189,0.485273,0.482284,0.479221,0.476084,0.472873,0.469588,0.466226,0.462789,0.459277,0.455688,0.452024,0.448284,0.444467,0.440573,0.436601,0.432552,0.428426,0.424223,0.419943,0.415586,0.411152,0.406640,0.402049,0.397381,0.392636,0.387814,0.382914,0.377939,0.372886,0.367757,0.362552,0.357269,0.351910,0.346476,0.340967,0.335384,0.329727,0.323998,0.318195,0.312321,0.306377,0.300362,0.294279,0.288127,0.281908,0.275626,0.269281,0.262877,0.256415,0.249897,0.243329,0.236712,0.230052,0.223353,0.216620,0.209861,0.203082,0.196293,0.189503,0.182725,0.175971,0.169257,0.162603,0.156029,0.149561,0.143228,0.137064,0.131109,0.125405,0.120005,0.114965,0.110347,0.106217,0.102646,0.099702,0.097452,0.095953,0.095250,0.095374,0.096335,0.098125,0.100717,0.104071,0.108131,0.112838,0.118128,0.123941,0.130215,0.136897,0.143936};
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float xmin,xmax,ymin,ymax,zmin,zmax=1.0;
  int i,j,k,flag=0,isel=0,sn;
  int redraw=1,mode=0,posn=0,click=0,graves=0,grid=0;
  float dt,zzmax,s1,s2,z,za,sigma,zs,zm;
  int ix=0,iy=0,isub=0;
  int i0,j0,i1,j1,jmax;
  float width=1500;
  float x=0.0,y=0.0,x0=0.0,y0=0.0,yfit;
  char c;
  char path[128],xlabel[128],ylabel[64],filename[32],tlefile[128];
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
  int cmap=2;
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

  // Set selection
  isel=0;
  sel.n=0;
  sel.w=100.0;
  sel.sigma=5.0;
  
  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:hc:C:gm:o:S:W:"))!=-1) {
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

      case 'W':
	sel.w=atof(optarg);
	break;

      case 'S':
	sel.sigma=atof(optarg);
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
  //  cpgpap(12.5,0.55);
  cpgask(0);

  // Default limits
  xmin=0.0;
  xmax=(float) s.nsub;
  ymin=0.0;
  ymax=(float) s.nchan;
  zmin=0.0;

  // Set trace
  tf.n=0;
  
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
      
      if (cmap==3) {
	cpggray(s.z,s.nsub,s.nchan,1,s.nsub,1,s.nchan,zmax,zmin,tr);
      } else {
	if (cmap==0)
	  cpgctab(cool_l,cool_r,cool_g,cool_b,9,1.0,0.5);
	else if (cmap==1)
	  cpgctab(heat_l,heat_r,heat_g,heat_b,9,1.0,0.5);
	else if (cmap==2)
	  cpgctab(viridis_l,viridis_r,viridis_g,viridis_b,256,1.0,0.5);
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
      
      // Human readable frequency axis
      fcen=0.5*(fmax+fmin);
      cpgswin(xmin,xmax,fmin-fcen,fmax-fcen);
      if (foverlay==1) 
	plot_traces(t,nsat,fcen);

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

    // Help
    if (c=='h') {
      printf("q        Quit\n");
      printf("k        Plot grid (not sure what it does)\n");
      printf("t        Locate traces\n");
      printf("s        Select track points\n");
      printf("u        Undo track point selection\n");
      printf("S        Set track selection width and sigma\n");
      printf("f        Fits selected track points (generates out.dat)\n");
      printf("F        Fits selected wideband track points (generates out.dat)\n");
      printf("g        Filter points above a certain sigma limit (generates filter.dat)\n");
      printf("G        Find peaks in current window (generates peakfind.dat)\n");
      printf("i        Identify selected track points\n");
      printf("I        Manually identify selected track points\n");
      printf("v        Decrease dynamic range\n");
      printf("b        Increase dynamic range\n");
      printf("Z        Provide manual dynamic range limits\n");
      printf("D/mid    Mark point (appends to mark.dat)\n");
      printf("w        Fit wideband point (appends to mark.dat)\n");
      printf("c        Center view\n");
      printf("C        Toggle through color maps\n");
      printf("p/right  Toggle overlays\n");
      printf("+        Zoom\n");
      printf("-/x      Unzoom\n");
      printf("R        Recompute traces\n");
      printf("r        Reset view\n");
      printf("z        Select part to zoom into\n");
      printf("TAB      Pan through view\n");
      
      continue;
    }
    
    
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
	if (graves==0)
	  locate_trace(s,t[i],site_id,sel.sigma,sel.w,sel.w,graves);
	else
	  locate_trace(s,t[i],site_id,sel.sigma,10.0,sel.w,graves);
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
      filter(s,site_id,sel.sigma,graves);

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

    // Fit
    if (c=='F') {
      tf=fit_gaussian_trace(s,sel,site_id,graves);
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

    // Set selection
    if (c=='S') {
      printf("Current selection width: %g pixels\n",sel.w);
      printf("Current sigma: %g\n",sel.sigma);
      printf("Provide selection width (pixels), and sigma: ");
      status=scanf("%f %f",&sel.w,&sel.sigma);
      printf("New selection width: %g pixels\n",sel.w);
      printf("New sigma: %g\n",sel.sigma);
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

    if (c=='Z') {
      printf("zmin,zmax\n");
      status=scanf("%f %f",&zmin,&zmax);
      printf("%f %f\n",zmin,zmax);
      redraw=1;
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
    if (c=='D' || c=='a') {
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

    // Fit single wideband point
    if (c=='w') {
      file=fopen("mark.dat","a");

      // Fit point
      yfit=fit_gaussian_point(s,x,y,sel,site_id,graves);
      cpgpt1(x,yfit,17);
      i=(int) floor(x);
      f=s.freq-0.5*s.samp_rate+(double) yfit*s.samp_rate/(double) s.nchan;
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
      if (cmap>3)
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
  free(s.length);
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
  free(gd.x);
  free(gd.y);
  free(gd.ym);

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
  printf("-o <offset>  Frequency offset to apply (Hz) [0]\n");
  printf("-W <width>   Track selection width (in pixels) [100]\n");
  printf("-S <sigma>   Track selection significance [5]\n");
  printf("-C <site>    Site ID\n");
  printf("-c <catalog> TLE catalog\n");
  printf("-g           GRAVES data\n");
  printf("-h           This help\n");

  return;
}


void plot_traces(struct trace *t,int nsat,float fcen)
{
  int i,j,flag,textflag;
  char text[8];

  // Loop over objects
  for (i=0;i<nsat;i++) {
    // Select color
    if (t[i].classfd==1)
      cpgsci(8);
    else
      cpgsci(3);
    
    sprintf(text," %d",t[i].satno);

    // Plot label at start of trace
    if (t[i].za[0]<=90.0)
	cpgtext(0.0,(float) t[i].freq[0]-fcen,text);

    // Loop over trace
    for (j=0,flag=0,textflag=0;j<t[i].n;j++) {
      // Plot label for rising sources
      if (j>0 && t[i].za[j-1]>90.0 && t[i].za[j]<=90.0)
	cpgtext((float) j,(float) t[i].freq[j]-fcen,text);

      // Plot line
      if (flag==0) {
	cpgmove((float) j,t[i].freq[j]-fcen);
	flag=1;
      } else {
	cpgdraw((float) j,t[i].freq[j]-fcen);
      }

      // Below horizon
      if (t[i].za[j]>90.0)
	flag=0;
      else
	flag=1;
    }	  
  }
  cpgsci(1);

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

      // Remove 
      s1-=zm;
      s2-=zm*zm;
      sn--;
      za=s1/(float) sn;
      zs=sqrt(s2/(float) sn-za*za);
      sigma=(zm-za)/zs;

      // Store
      if (sigma>sel.sigma && s.mjd[i]>1.0) {
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

// Fit gaussian trace
struct trace fit_gaussian_trace(struct spectrogram s,struct select sel,int site_id,int graves)
{
  int i,j,k,l,sn;
  int i0,i1,j0,j1,jmax;
  double f,chi2;
  double a[4],da[4];
  float x,y,s1,s2,z,za,zs,zm,sigma,rms;
  struct trace t;
  FILE *file;

  // Set up trace
  t.satno=99999;
  t.n=(int) ceil(sel.x[sel.n-1]-sel.x[0]);
  t.mjd=(double *) malloc(sizeof(double)*t.n);
  t.freq=(double *) malloc(sizeof(double)*t.n);
  t.za=(float *) malloc(sizeof(float)*t.n);

  // Set up data
  gd.n=(int) 2 * sel.w;
  gd.x=(double *) malloc(sizeof(double) * gd.n);
  gd.y=(double *) malloc(sizeof(double) * gd.n);
  gd.ym=(double *) malloc(sizeof(double) * gd.n);
  
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

      // Loop over points
      for (j=j0,l=0;j<j1;j++,l++) {
	gd.x[l] = (double) j;
	gd.y[l] = (double) s.z[i+s.nsub*j];
      }
      gd.n = l;

      // Parameter guess
      a[0]=0.5*(j0+j1);
      da[0]=10.0;
      a[1]=10.0;
      da[1]=1.0;
      a[2]=1.0;
      da[2]=1.0;
      a[3]=1.0;
      da[3]=1.0;

      // Run fit
      versafit(gd.n,4,a,da,chisq_gaussian,0.0,1e-3,"n");
      
      // Find maximum and significance
      for (j=0,s1=0.0;j<gd.n;j++) {
	s1+=pow(gd.y[i]-gd.ym[i],2);
      }
      zs=sqrt(s1/(float) gd.n);
      sigma=a[2]/zs;
      
      // Store
      if (sigma>sel.sigma && s.mjd[i]>1.0) {
	f=s.freq-0.5*s.samp_rate+(double) a[0]*s.samp_rate/(double) s.nchan;
	if (graves==0)
	  fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sigma,site_id);
	else
	  fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,sigma,site_id);
	cpgpt1((float) i,(float) a[0], 17);
	cpgmove((float) i,(float) a[0]-a[1]);
	cpgdraw((float) i,(float) a[0]+a[1]);
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

// Fit trace
struct trace locate_trace(struct spectrogram s,struct trace t,int site_id,float sigmamin,float wmin,float wmax,int graves)
{
  int i,j,k,l,sn;
  int i0,i1,j0,j1,jmax;
  double f,fmin;
  float x,y,s1,s2,z,za,zs,zm,sigma;
  FILE *file;
  char filename[64];
  int file_open=0;
  

  fmin=(s.freq-0.5*s.samp_rate)*1e-6;

  // Loop over trace
  for (i=0;i<t.n;i++) {
    // Skip when satellite is below the horizon
    if (t.za[i]>90.0)
      continue;
    
    // Compute position
    y=(t.freq[i]-fmin)*s.nchan/(s.samp_rate*1e-6);
    j0=(int) floor(y-wmax);
    j1=(int) floor(y+wmax);

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
    if (sigma>sigmamin && s.mjd[i]>1.0 && fabs(y-jmax)<wmin) {
      // Open file
      if (file_open==0) {
	// Open file
	sprintf(filename,"track_%05d_%08.3f.dat",t.satno,t.freq0);
	file=fopen(filename,"a");
	file_open=1;
      }
      f=s.freq-0.5*s.samp_rate+(double) jmax*s.samp_rate/(double) s.nchan;
      if (graves==0)
	fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sigma,site_id);
      else
	fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,sigma,site_id);
      cpgpt1((float) i,(float) jmax,17);
    }
  }

  // Close file
  if (file_open==1)
    fclose(file);

  return t;
}

// Filter data
void filter(struct spectrogram s,int site_id,float sigma,int graves)
{
  int i,j,k,l,jmax,zmax;
  float s1,s2,avg,std,dz;
  FILE *file;
  double f;
  int *mask;

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

    // Mark points
    for (j=0;j<s.nchan;j++) {
      if (mask[j]==1) {
	f=s.freq-0.5*s.samp_rate+(double) j*s.samp_rate/(double) s.nchan;
	if (s.mjd[i]>1.0) {
	  if (graves==0)
	    fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id);
	  else
	    fprintf(file,"%lf %lf %f %d 9999 %d %d\n",s.mjd[i],f,s.z[i+s.nsub*j],site_id,i,j);
	}
	cpgpt1((float) i+0.5,(float) j+0.5,17);
      }
    }
  }
  fclose(file);

  free(mask);

  return;
}

// Peakfinding algorithm
void peakfind(struct spectrogram s,int site_id,int i0,int i1,int j0,int j1)
{
  int i,j,k,l,m=21,n;
  float *w,*y,*sy,*a,*b,*c,d[3],dx[3],dw=1.0,x0,c0=-0.0007;
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

// Gaussian point
double chisq_gaussian(double a[])
{
  int i;
  double arg,amp,chisq;

  // Compute chisq
  for (i=0,chisq=0.0;i<gd.n;i++) {
    // Compute gaussian
    arg=-0.5*pow((gd.x[i]-a[0])/a[1],2);
    amp=a[2];
    gd.ym[i]=amp*exp(arg)+a[3];

    // Sum to chi-squared
    chisq+=pow(gd.y[i]-gd.ym[i],2);
  }
    
  return chisq;
}

// Fit gaussian point
float fit_gaussian_point(struct spectrogram s,float x,float y,struct select sel,int site_id,int graves)
{
  int i,j,k,l,sn;
  int i0,i1,j0,j1,jmax;
  double f,chi2;
  double a[4],da[4];

  // Set up data
  gd.n=(int) 2 * sel.w;
  gd.x=(double *) malloc(sizeof(double) * gd.n);
  gd.y=(double *) malloc(sizeof(double) * gd.n);
  gd.ym=(double *) malloc(sizeof(double) * gd.n);
  
  // Set point
  i=(int) floor(x);
  j0=(int) floor(y-sel.w);
  j1=(int) floor(y+sel.w);

  // Keep in range
  if (j0<0)
    j0=0;
  if (j1>=s.nchan)
    j1=s.nchan;
  
  // Loop over points
  for (j=j0,l=0;j<j1;j++,l++) {
    gd.x[l] = (double) j;
    gd.y[l] = (double) s.z[i+s.nsub*j];
  }
  gd.n = l;

  // Parameter guess
  a[0]=0.5*(j0+j1);
  da[0]=10.0;
  a[1]=10.0;
  da[1]=1.0;
  a[2]=1.0;
  da[2]=1.0;
  a[3]=1.0;
  da[3]=1.0;

  // Run fit
  versafit(gd.n,4,a,da,chisq_gaussian,0.0,1e-3,"n");
      
  return a[0];
}
