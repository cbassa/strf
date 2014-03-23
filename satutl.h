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
int  read_twoline(FILE *fp, long satno, orbit_t *orb);
void *vector(size_t num, size_t size);
void print_orb(orbit_t *orb);

/** aries.c **/
double gha_aries(double jd);

/** ferror.c **/
void fatal_error(const char *format, ...);

#ifdef __cplusplus
}
#endif

#endif /* _SATUTL_H */
