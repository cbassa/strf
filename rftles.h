#ifndef _RFTLES_H
#define _RFTLES_H

#include "sgdp4h.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tles {
    long number_of_elements;
    orbit_t *orbits;
} tles_t;

tles_t load_tles(char *tlefile);
void free_tles(tles_t *tles);
orbit_t *get_orbit_by_index(tles_t *tles, long index);
orbit_t *get_orbit_by_catalog_id(tles_t *tles, long satno);

#ifdef __cplusplus
}
#endif

#endif /* _RFTLES_H */
