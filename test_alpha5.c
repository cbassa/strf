#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rftles.h"

#define LIM 80

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

int main(int argc,char *argv[])
{
  tle_t *tle;
  char tlefile[]="alpha5_test.txt";
  
  // Load TLEs
  tle_array_t *tle_array = load_tles(tlefile);

  if (tle_array->number_of_elements == 0) {
    fprintf(stderr,"TLE file %s not found or empty\n", tlefile);
    return 0;
  }

  // Loop over all TLEs
  for (long elem = 0; elem < tle_array->number_of_elements; elem++) {
    // Get TLE
    tle = get_tle_by_index(tle_array, elem);

    print_orb(&tle->orbit);
    printf("\n");
  }

  // Find specific TLE
  tle = get_tle_by_catalog_id(tle_array, 5);
  print_orb(&tle->orbit);
  printf("\n");

  // Find specific TLE
  tle = get_tle_by_catalog_id(tle_array, 115544);
  print_orb(&tle->orbit);
  printf("\n");
  

  return 0;
}
