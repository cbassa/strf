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
struct point {
  double mjd;
  double freq;
};

void dec2sex(double x,char *s,int f,int len);
void time_axis(double *mjd,int n,float xmin,float xmax,float ymin,float ymax);
void usage(void);
void plot_traces(struct trace *t,int nsat,float foff);
void filter(struct spectrogram s,int site_id,float sigma,char *filename,int graves);

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
  float xmin,xmax,ymin,ymax,zmin,zmax=8.0;
  int i,j,k;
  float dt,zzmax,s1,s2;
  int ix=0,iy=0,isub=0;
  int i0,j0,i1,j1,jmax;
  float width=1500,sigma=5.0,foff=0.0;
  float x,y,x0,y0;
  char c;
  char path[128],xlabel[64],ylabel[64],filename[32],tlefile[128],pngfile[128],datfile[128];
  int sec,lsec,ssec;
  char stime[16];
  double fmin,fmax,fcen,f;
  FILE *file;
  int arg=0,nsub=1800,nbin=1;
  double f0=0.0,df0=0.0,dy=2500;
  int foverlay=1;
  struct trace *t,tf;
  int nsat,satno;
  struct select sel;
  char *env;
  int site_id=0,cmap=2,graves=0,create_dat=1,fname_flag=1;

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
    while ((arg=getopt(argc,argv,"p:f:w:s:l:b:z:hc:C:m:gS:qo:O:"))!=-1) {
      switch (arg) {
	
      case 'p':
	strcpy(path,optarg);
	break;

      case 'o':
	strcpy(pngfile,optarg);
	fname_flag=0;
	break;

      case 'O':
	foff=atof(optarg);
	break;
	
      case 's':
	isub=atoi(optarg);
	break;
	
      case 'f':
	f0=(double) atof(optarg);
	break;

      case 'l':
	nsub=atoi(optarg);
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

      case 'g':
	graves=1;
	break;

      case 'S':
	sigma=atof(optarg);
	break;

      case 'q':
	create_dat=0;
	break;
	
      case 'C':
	site_id=atoi(optarg);
	break;
	
      case 'm':
	cmap=atoi(optarg);
	if (cmap>2)
	  cmap=0;
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

  // Read data
  s=read_spectrogram(path,isub,nsub,f0,df0,nbin,foff);
  if (s.mjd[0]<54000)
    return 0;

  // Output filename
  if (fname_flag==1)
    sprintf(pngfile,"%.19s_%08.3f.png/png",s.nfd0,s.freq*1e-6);

  if (create_dat==1)
    sprintf(datfile,"%.19s_%08.3f.dat",s.nfd0,s.freq*1e-6);
  
  printf("Read spectrogram\n%d channels, %d subints\nFrequency: %g MHz\nBandwidth: %g MHz\n",s.nchan,s.nsub,s.freq*1e-6,s.samp_rate*1e-6);

  // Compute traces
  t=compute_trace(tlefile,s.mjd,s.nsub,site_id,s.freq*1e-6,s.samp_rate*1e-6,&nsat,graves);
  printf("Traces for %d objects for location %d\n",nsat,site_id);

  cpgopen(pngfile);
  //    cpgctab(cool_l,cool_r,cool_g,cool_b,9,1.0,0.5);
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
    
  cpgimag(s.z,s.nsub,s.nchan,1,s.nsub,1,s.nchan,zmin,zmax,tr);
    
  // Pixel axis
  cpgbox("CTSM1",0.,0,"CTSM1",0.,0);
  
  // Time axis
  cpgbox("B",0.,0,"",0.,0);

  cpgswin(xmin,xmax,ymin,ymax);

  // Filter points
  if (create_dat)
    filter(s,site_id,sigma,datfile,graves);
    
  time_axis(s.mjd,s.nsub,xmin,xmax,ymin,ymax);
  
  // Freq axis
  fmin=s.freq-0.5*s.samp_rate+ymin*s.samp_rate/(float) s.nchan;
  fmax=s.freq-0.5*s.samp_rate+ymax*s.samp_rate/(float) s.nchan;
  fmin*=1e-6;
  fmax*=1e-6;
  
  // Plot traces
  cpgswin(xmin,xmax,fmin,fmax);
  cpgsch(0.8);
  plot_traces(t,nsat,0.0);
  cpgsch(0.8);
  
  // Human readable frequency axis
  fcen=0.5*(fmax+fmin);
  fcen=floor(1000*fcen)/1000.0;
  if (fabs(foff)<1e-3)
    sprintf(ylabel,"Frequency - %.3f MHz",fcen);
  else
    sprintf(ylabel,"Frequency - %.3f MHz (%.3f Hz)",fcen, foff);
  fmin-=fcen;
  fmax-=fcen;
  cpgswin(xmin,xmax,fmin,fmax);
  cpgbox("",0.,0,"BTSN",0.,0);
    
  sprintf(xlabel,"UT Date: %.10s",s.nfd0);
  cpglab(xlabel,ylabel," ");
  
  
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
  //  lsec=60;
  //  ssec=10;

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
  printf("-O <offset>  Frequency offset to apply (Hz) [0]\n");
  printf("-h           This help\n");

  return;
}

void plot_traces(struct trace *t,int nsat,float foff)
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
	cpgtext(0.0,(float) t[i].freq[0],text);

    // Loop over trace
    for (j=0,flag=0,textflag=0;j<t[i].n;j++) {
      // Plot label for rising sources
      if (j>0 && t[i].za[j-1]>90.0 && t[i].za[j]<=90.0)
	cpgtext((float) j,(float) (t[i].freq[j]+foff),text);

      // Plot line
      if (flag==0) {
	cpgmove((float) j,(float) (t[i].freq[j]+foff));
	flag=1;
      } else {
	cpgdraw((float) j,(float) (t[i].freq[j]+foff));
      }

      // Below horizon
      if (t[i].za[j]>90.0)
	flag=0;
      else
	flag=1;
    }
    cpgsci(1);
  }

  return;
}

// Filter points
void filter(struct spectrogram s,int site_id,float sigma,char *filename,int graves)
{
  int i,j,k,l;
  float s1,s2,avg,std,dz;
  FILE *file;
  double f;
  int *mask;
  float *sig;

  mask=(int *) malloc(sizeof(int)*s.nchan);
  sig=(float *) malloc(sizeof(float)*s.nchan);
  
  // Open file
  file=fopen(filename,"w");

  // Loop over subints
  for (i=0;i<s.nsub;i++) {
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
      sig[j]=(s.z[i+s.nsub*j]-avg)/std;
      if (sig[j]>sigma) 
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
	    fprintf(file,"%lf %lf %f %d\n",s.mjd[i],f,sig[j],site_id);
	  else
	    fprintf(file,"%lf %lf %f %d 9999\n",s.mjd[i],f,sig[j],site_id);
	}
	cpgpt1((float) i+0.5,(float) j+0.5,17);
      }
    }
  } 

  fclose(file);

  free(mask);
  free(sig);

  return;
}

