#include "rftles.h"

#include "satutl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

tle_array_t load_tles(char *tlefile) {
    tle_array_t tle_array;

    tle_array.tles = NULL;
    tle_array.number_of_elements = 0;

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

    size_t linesize = 256;
    char * line = malloc(linesize);
    ssize_t read;

    // Count number of entries
    long num_elements = 0;

    while ((read = getline(&line, &linesize, file)) != -1) {
        if (read > 0 && strncmp(line, "1 ", 2) == 0) {
            num_elements++;
        }
    }

    // Don't allocate anything if no entry found in file
    if (num_elements == 0) {
      fclose(file);
      free(line);

      return tle_array;
    }

    tle_array.tles = (tle_t *)calloc(num_elements, sizeof(tle_t));

    // Rewind and parse file
    rewind(file);

    while (read_twoline(file, 0, &(tle_array.tles[tle_array.number_of_elements].orbit)) == 0) {
        tle_array.number_of_elements++;
    }

    free(line);
    fclose(file);

    printf("Loaded %ld orbits\n", tle_array.number_of_elements);

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
        tle_array->number_of_elements = 0;
    }
}

tle_t *get_orbit_by_index(tle_array_t *tle_array, long index) {
    if (tle_array && (index < tle_array->number_of_elements)) {
        return &(tle_array->tles[index]);
    }

    return NULL;
}

tle_t *get_orbit_by_catalog_id(tle_array_t *tle_array, long satno) {
    if (tle_array) {
        for (long i = 0; i < tle_array->number_of_elements; i++) {
            if (tle_array->tles[i].orbit.satno == satno) {
                return &(tle_array->tles[i]);
            }
        }
    }

    return NULL;
}
