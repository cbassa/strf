#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc,char *argv[])
{
  int i,firstfile=1,status;
  FILE *file;
  char header[256],filename[128],nfd[32],fileroot[128];
  double freq,samp_rate;
  float length;
  int nchan;

  if (strchr(argv[1],'_')!=NULL) {
    strncpy(fileroot,argv[1],strlen(argv[1])-11);
    fileroot[strlen(argv[1])-11]='\0';
  } else {
    strcpy(fileroot,argv[1]);
  }

  // Open first file
  for (i=0;;i++) {
    sprintf(filename,"%s_%06d.bin",fileroot,i);
    file=fopen(filename,"r");
    
    // Break if file does not exist
    if (file==NULL) 
      break;

    // Read header
    if (firstfile==1) {
      // Read header
      status=fread(header,sizeof(char),256,file);
      status=sscanf(header,"HEADER\nUTC_START    %s\nFREQ         %lf Hz\nBW           %lf Hz\nLENGTH       %f s\nNCHAN        %d\n",nfd,&freq,&samp_rate,&length,&nchan);
      firstfile=0;
    }

    fclose(file);
  }

  printf("%s %8.3lf %8.3lf %d %d\n",argv[1],freq*1e-6,samp_rate*1e-6,nchan,i);

  return 0;
}
