#include "rftles.h"

#include "elset.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

tle_array_t *load_tles(char *tlefile) {
    tle_array_t *tle_array;

    tle_array = (tle_array_t *)malloc(sizeof(tle_array_t));

    if (tle_array == NULL) {
      return NULL;
    }

    tle_array->tles = NULL;
    tle_array->number_of_elements = 0;

    char filename[1024];

    if (tlefile) {
        strncpy(filename, tlefile, sizeof(filename));
    } else {
        char * env = getenv("ST_TLEDIR");

        if (env == NULL || strlen(env) == 0) {
            env=".";
        }

        sprintf(filename, "%s/bulk.tle", env);
    }

    FILE * file = fopen(filename, "r");

    if (file == NULL) {
        fprintf(stderr, "TLE file %s not found\n", filename);

        return tle_array;
    }

    // Read element sets in any format read_elset() knows (TLE, 3LE,
    // OMM as JSON, XML, KVN or CSV), growing the array as needed
    long capacity = 0;
    orbit_t orbit;
    char satname[ST_SIZE];

    memset(&orbit, 0, sizeof(orbit));
    while (read_elset(file, 0, &orbit, satname) == 0) {
        if (tle_array->number_of_elements == capacity) {
            long new_capacity = (capacity == 0) ? 256 : 2 * capacity;
            tle_t *tles = (tle_t *)realloc(tle_array->tles, new_capacity * sizeof(tle_t));

            if (tles == NULL) {
                fclose(file);
                free_tles(tle_array);

                return NULL;
            }
            tle_array->tles = tles;
            capacity = new_capacity;
        }

        tle_t *tle = &(tle_array->tles[tle_array->number_of_elements]);

        tle->orbit = orbit;
        if (satname[0] != '\0') {
            int satname_len = strlen(satname) + 1;
            tle->name = malloc(satname_len);
            strncpy(tle->name, satname, satname_len);
        } else {
            tle->name = NULL;
        }

        tle_array->number_of_elements++;
        memset(&orbit, 0, sizeof(orbit));
    }

    fclose(file);

    printf("Loaded %ld orbits\n", tle_array->number_of_elements);

    return tle_array;
}

void free_tles(tle_array_t *tle_array) {
    if (tle_array) {
        for (long i = 0; i < tle_array->number_of_elements; i++) {
            if (tle_array->tles[i].name != NULL) {
                free(tle_array->tles[i].name);
            }
        }

        free(tle_array->tles);
        free(tle_array);
    }
}

tle_t *get_tle_by_index(tle_array_t *tle_array, long index) {
    if (tle_array && (index < tle_array->number_of_elements)) {
        return &(tle_array->tles[index]);
    }

    return NULL;
}

tle_t *get_tle_by_catalog_id(tle_array_t *tle_array, long satno) {
    if (tle_array) {
        for (long i = 0; i < tle_array->number_of_elements; i++) {
            if (tle_array->tles[i].orbit.satno == satno) {
                return &(tle_array->tles[i]);
            }
        }
    }

    return NULL;
}
