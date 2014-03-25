#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>

void usage(void)
{
  printf("rffft: FFT RF observations\n\n");
  printf("-i <file>       Input file (can be fifo)\n");
  printf("-p <prefix>     Output prefix\n");
  printf("-f <frequency>  Center frequency (Hz)\n");
  printf("-s <samprate>   Sample rate (Hz)\n");
  printf("-c <chansize>   Channel size [100Hz]\n");
  printf("-t <tint>       Integration time [1s]\n");
  printf("-n <nsub>       Number of integrations per file [60]\n");
  printf("-h              This help\n");

  return;
}

int main(int argc,char *argv[])
{
  int i,j,k,l,m,nchan,nint=1,arg=0,nbytes,nsub=60,flag;
  fftwf_complex *c,*d;
  fftwf_plan fft;
  FILE *infile,*outfile;
  char infname[128],outfname[128],path[64]="/data/record",prefix[32]="";
  float *buf,tint=1.0;
  float *z,length,fchan=100.0;
  double freq,samp_rate;
  struct timeval start,end;
  char tbuf[30],nfd[32],header[256]="";

  // Read arguments
  if (argc>1) {
    while ((arg=getopt(argc,argv,"i:f:s:c:t:p:n:h"))!=-1) {
      switch(arg) {
	
      case 'i':
	strcpy(infname,optarg);
	break;
	
      case 'p':
	strcpy(path,optarg);
	break;
	
      case 'f':
	freq=(double) atof(optarg);
	break;
	
      case 's':
	samp_rate=(double) atof(optarg);
	break;
	
      case 'c':
	fchan=atof(optarg);
	break;
	
      case 'n':
	nsub=atoi(optarg);
	break;
	
      case 't':
	tint=atof(optarg);
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

  // Number of channels
  nchan=(int) (samp_rate/fchan);

  // Number of integrations
  nint=(int) (tint*(float) samp_rate/(float) nchan);

  // Dump statistics
  printf("Filename: %s\n",infname);
  printf("Bandwidth: %f MHz\n",samp_rate*1e-6);
  printf("Sampling time: %f us\n",1e6/samp_rate);
  printf("Number of channels: %d\n",nchan);
  printf("Channel size: %f Hz\n",samp_rate/(float) nchan);
  printf("Integration time: %f s\n",tint);
  printf("Number of averaged spectra: %d\n",nint);
  printf("Number of subints per file: %d\n",nsub);

  // Allocate
  c=fftwf_malloc(sizeof(fftwf_complex)*nchan);
  d=fftwf_malloc(sizeof(fftwf_complex)*nchan);
  buf=(float *) malloc(sizeof(float)*2*nchan);
  z=(float *) malloc(sizeof(float)*nchan);

  // Plan
  fft=fftwf_plan_dft_1d(nchan,c,d,FFTW_FORWARD,FFTW_ESTIMATE);

  // Create prefix
  gettimeofday(&start,0);
  strftime(prefix,30,"%Y-%m-%dT%T",gmtime(&start.tv_sec));

  // Open file
  infile=fopen(infname,"r");

  // Forever loop
  for (m=0;;m++) {
    // File name
    sprintf(outfname,"%s/%s_%06d.bin",path,prefix,m);
    outfile=fopen(outfname,"w");

    // Loop over subints to dump
    for (k=0;k<nsub;k++) {
      // Initialize
      for (i=0;i<nchan;i++) 
	z[i]=0.0;
      
      // Log start time
      gettimeofday(&start,0);
      
      // Integrate
      for (j=0;j<nint;j++) {
	nbytes=fread(buf,sizeof(float),2*nchan,infile);
	if (nbytes==0)
	  break;
	
	// Unpack
	for (i=0;i<nchan;i++) {
	  c[i][0]=buf[2*i];
	  c[i][1]=buf[2*i+1];
	}
	
	// Execute
	fftwf_execute(fft);
	
	// Add
	for (i=0;i<nchan;i++) {
	  if (i<nchan/2)
	    l=i+nchan/2;
	  else
	    l=i-nchan/2;
	  
	  z[l]+=d[i][0]*d[i][0]+d[i][1]*d[i][1];
	}
      }

      // Log end time
      gettimeofday(&end,0);
      
      // Scale
      for (i=0;i<nchan;i++) 
	z[i]/=(float) nchan;
      
      // Time stats
      length=(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)*1e-6;
      
      // Format start time
      strftime(tbuf,30,"%Y-%m-%dT%T",gmtime(&start.tv_sec));
      sprintf(nfd,"%s.%03ld",tbuf,start.tv_usec/1000);

      // Header
      sprintf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\nNSUB         %d\nEND\n",nfd,freq,samp_rate,length,nchan,nsub);

      printf("%s %s %f %d\n",outfname,nfd,length,j);
      
      // Dump file
      fwrite(header,sizeof(char),256,outfile);
      fwrite(z,sizeof(float),nchan,outfile);
      
      // Break;
      if (nbytes==0)
	break;
    }

    // Break;
    if (nbytes==0)
      break;
    
    // Close file
    fclose(outfile);
  }
  fclose(infile);

  // Destroy plan
  fftwf_destroy_plan(fft);

  // Deallocate
  free(buf);
  fftwf_free(c);
  fftwf_free(d);
  free(z);

  return 0;
}

