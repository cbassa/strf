/* > satutl.c
 *
 */


#include "sgdp4h.h"

#include <ctype.h>

static char *st_start(char *buf);
static long i_read(char *str, int start, int stop);
static double d_read(char *str, int start, int stop);

 /* ====================================================================
   Read a string from key board, remove CR/LF etc.
   ==================================================================== */

void read_kb(char *buf)
{
int ii;

    fgets(buf, ST_SIZE-1, stdin);

    /* Remove the CR/LF etc. */
    for(ii = 0; ii < ST_SIZE; ii++)
        {
        if(buf[ii] == '\r' || buf[ii] == '\n')
            {
            buf[ii] = '\0';
            break;
            }
        }
}

/* ====================================================================
   Read orbit parameters for "satno" in file "filename", return -1 if
   failed to find the corresponding data. Call with satno = 0 to get the
   next elements of whatever sort.
   ==================================================================== */

int read_twoline(FILE *fp, long search_satno, orbit_t *orb)
{
  static char search[ST_SIZE];
  static char line1[ST_SIZE];
  static char line2[ST_SIZE];
  char *st1, *st2;
  int found;
  double bm, bx;
  
  st1 = line1;
  st2 = line2;
  
  do {
    if(fgets(line1, ST_SIZE-1, fp) == NULL) return -1;
    st1 = st_start(line1);
  } while(st1[0] != '1');
  
  if(search_satno > 0)
    {
      found = 0;
    }
  else
    {
      found = 1;
      search_satno = atol(st1+2);
    }
  
  sprintf(search, "1 %05ld", search_satno);
  
  do {
    st1 = st_start(line1);
    if(strncmp(st1, search, 7) == 0)
      {
	found = 1;
	break;
      }
  } while(fgets(line1, ST_SIZE-1, fp) != NULL);
  
  
  sprintf(search, "2 %05ld", search_satno);
  
  if(found)
    {
      fgets(line2, ST_SIZE-1, fp);
      st2 = st_start(line2);
    }
  
  if(!found || strncmp(st2, search, 7) != 0)
    {
      return -1;
    }
  
  orb->ep_year = (int)i_read(st1, 19, 20);
  
  if(orb->ep_year < 57) orb->ep_year += 2000;
  else orb->ep_year += 1900;
  
  orb->ep_day =       d_read(st1, 21, 32);
  
  orb->ndot2 = d_read(st1, 34, 43);
  bm = d_read(st1, 45, 50) * 1.0e-5;
  bx = d_read(st1, 51, 52);
  orb->nddot6 = bm * pow(10.0, bx);
  bm = d_read(st1, 54, 59) * 1.0e-5;
  bx = d_read(st1, 60, 61);
  orb->bstar = bm * pow(10.0, bx);

  orb->eqinc = RAD(d_read(st2,  9, 16));
  orb->ascn = RAD(d_read(st2, 18, 25));
  orb->ecc  =     d_read(st2, 27, 33) * 1.0e-7;
  orb->argp = RAD(d_read(st2, 35, 42));
  orb->mnan = RAD(d_read(st2, 44, 51));
  orb->rev  =     d_read(st2, 53, 63);
  orb->norb =     i_read(st2, 64, 68);
  
  orb->satno = search_satno;

  sscanf(st1+9,"%s",orb->desig);

  return 0;
}

/* ==================================================================
   Locate the first non-white space character, return location.
   ================================================================== */

static char *st_start(char *buf)
{
    if(buf == NULL) return buf;

    while(*buf != '\0' && isspace(*buf)) buf++;

return buf;
}

/* ==================================================================
   Mimick the FORTRAN formatted read (assumes array starts at 1), copy
   characters to buffer then convert.
   ================================================================== */

static long i_read(char *str, int start, int stop)
{
long itmp=0;
char *buf, *tmp;
int ii;

    start--;    /* 'C' arrays start at 0 */
    stop--;

    tmp = buf = (char *)vector(stop-start+2, sizeof(char));

    for(ii = start; ii <= stop; ii++)
        {
        *tmp++ = str[ii];   /* Copy the characters. */
        }
    *tmp = '\0';            /* NUL terminate */

    itmp = atol(buf);       /* Convert to long integer. */
    free(buf);

return itmp;
}

/* ==================================================================
   Mimick the FORTRAN formatted read (assumes array starts at 1), copy
   characters to buffer then convert.
   ================================================================== */

static double d_read(char *str, int start, int stop)
{
double dtmp=0;
char *buf, *tmp;
int ii;

    start--;
    stop--;

    tmp = buf = (char *)vector(stop-start+2, sizeof(char));

    for(ii = start; ii <= stop; ii++)
        {
        *tmp++ = str[ii];   /* Copy the characters. */
        }
    *tmp = '\0';            /* NUL terminate */

    dtmp = atof(buf);       /* Convert to long integer. */
    free(buf);

return dtmp;
}

/* ==================================================================
   Allocate and check an all-zero array of memory (storage vector).
   ================================================================== */

void *vector(size_t num, size_t size)
{
void *ptr;

    ptr = calloc(num, size);
    if(ptr == NULL)
        {
        fatal_error("vector: Allocation failed %u * %u", num, size);
        }

return ptr;
}

/* ==================================================================
   Print out orbital parameters.
   ================================================================== */

void print_orb(orbit_t *orb)
{
  printf("# Satellite ID = %ld\n", (long)orb->satno);
  printf("# Satellite designation = %s\n",orb->desig);
  printf("# Epoch year = %d day = %.8f\n", orb->ep_year, orb->ep_day);
  printf("# Eccentricity = %.7f\n", orb->ecc);
  printf("# Equatorial inclination = %.4f deg\n", DEG(orb->eqinc));
  printf("# Argument of perigee = %.4f deg\n", DEG(orb->argp));
  printf("# Mean anomaly = %.4f deg\n", DEG(orb->mnan));
  printf("# Right Ascension of Ascending Node = %.4f deg\n", DEG(orb->ascn));
  printf("# Mean Motion (number of rev/day) = %.8f\n", orb->rev);
  printf("# First derivative of mean motion = %e\n",orb->ndot2);
  printf("# Second derivative of mean motion = %e\n",orb->nddot6);
  printf("# BSTAR drag = %.4e\n", orb->bstar);
  printf("# Orbit number = %ld\n", orb->norb);
}

/* ====================================================================== */
