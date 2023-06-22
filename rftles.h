#ifndef _RFTLES_H
#define _RFTLES_H

#include "sgdp4h.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tle {
    orbit_t orbit;
    char *name;
} tle_t;

typedef struct tle_array {
    long number_of_elements;
    tle_t *tles;
} tle_array_t;

tle_array_t load_tles(char *tlefile);
void free_tles(tle_array_t *tle_array);
tle_t *get_tle_by_index(tle_array_t *tle_array, long index);
tle_t *get_tle_by_catalog_id(tle_array_t *tle_array, long satno);

#ifdef __cplusplus
}
#endif

#endif /* _RFTLES_H */
