/* > satutl.h
 *
 */

#ifndef _SATUTL_H
#define _SATUTL_H

#define ST_SIZE 256

#ifdef __cplusplus
extern "C" {
#endif

/** satutl.c **/
void read_kb(char *buf);
void conditional_copy_satname(char * satname, char * current_line);
int  read_twoline(FILE *fp, long satno, orbit_t *orb, char *satname);
void *vector(size_t num, size_t size);
void print_orb(orbit_t *orb);
int alpha5_to_number(const char *s);
void number_to_alpha5(int number, char *result);
void zero_pad(const char *s, char *result);
void strip_leading_spaces(const char *s, char *result);
void strip_trailing_spaces(char *s);
  
/** aries.c **/
double gha_aries(double jd);

/** ferror.c **/
void fatal_error(const char *format, ...);
  
#ifdef __cplusplus
}
#endif

#endif /* _SATUTL_H */
