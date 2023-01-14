#ifndef _RFSITES_H
#define _RFSITES_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct site {
  int id;
  double lat,lng;
  float alt;
  char observer[64];
} site_t;

site_t get_site(int site_id);

#ifdef __cplusplus
}
#endif

#endif /* _RFSITES_H */
