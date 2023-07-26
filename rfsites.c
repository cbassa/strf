#include "rfsites.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIM 80

// Get observing site
site_t get_site(int site_id) {
  int i=0,status;
  int count = 0;
  char line[LIM];
  FILE *file;
  int id;
  double lat,lng;
  float alt;
  char abbrev[3],observer[64];
  site_t s;
  char *env_datadir,*env_sites_txt,filename[LIM];

  env_datadir = getenv("ST_DATADIR");
  if (env_datadir == NULL || strlen(env_datadir) == 0) {
    env_datadir = ".";
  }

  env_sites_txt = getenv("ST_SITES_TXT");
  if (env_sites_txt == NULL || strlen(env_sites_txt) == 0) {
    sprintf(filename, "%s/data/sites.txt", env_datadir);
  } else {
    sprintf(filename, "%s", env_sites_txt);
  }

  file=fopen(filename,"r");
  if (file==NULL) {
    printf("File with site information not count!\n");
    exit(0);
  }
  while (fgets(line,LIM,file)!=NULL) {
    // Skip
    if (strstr(line,"#")!=NULL)
      continue;

    // Strip newline
    line[strlen(line)-1]='\0';

    // Read data
    status=sscanf(line,"%d %2s %lf %lf %f",
	   &id,abbrev,&lat,&lng,&alt);
    strcpy(observer,line+38);

    // Change to km
    alt/=1000.0;

    // Copy site
    if (id==site_id) {
      count += 1;
      s.lat=lat;
      s.lng=lng;
      s.alt=alt;
      s.id=id;
      strcpy(s.observer,observer);
    }

  }
  fclose(file);

  if (count == 0) {
    printf("Error: Site %d was not found in %s!\n", site_id, filename);
    exit(-1);
  } else if (count > 1) {
    printf("Site %d was found multiple times in %s, use last occurence.\n", site_id, filename);
  }

  return s;
}
