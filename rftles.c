#include "rftles.h"

#include "satutl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

tles_t load_tles(char *tlefile) {
    tles_t tles;

    tles.orbits = NULL;
    tles.number_of_elements = 0;

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

        return tles;
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

      return tles;
    }

    tles.orbits = (orbit_t *)calloc(num_elements, sizeof(orbit_t));

    // Rewind and parse file
    rewind(file);

    while (read_twoline(file, 0, &(tles.orbits[tles.number_of_elements])) == 0) {
        tles.number_of_elements++;
    }

    free(line);
    fclose(file);

    printf("Loaded %ld orbits\n", tles.number_of_elements);

    return tles;
}

void free_tles(tles_t *tles) {
    if (tles) {
        free(tles->orbits);
        tles->number_of_elements = 0;
    }
}

orbit_t *get_orbit_by_index(tles_t *tles, long index) {
    if (tles && (index < tles->number_of_elements)) {
        return &(tles->orbits[index]);
    }

    return NULL;
}

orbit_t *get_orbit_by_catalog_id(tles_t *tles, long satno) {
    if (tles) {
        for (long i = 0; i < tles->number_of_elements; i++) {
            if (tles->orbits[i].satno == satno) {
                return &(tles->orbits[i]);
            }
        }
    }

    return NULL;
}
