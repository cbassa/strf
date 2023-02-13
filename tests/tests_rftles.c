#include "tests_rftles.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmocka.h>

#include "../rftles.h"

// Helpers functions for setup/teardown
void *tle_load(char *filename) {
  tle_array_t *tle_array = malloc(sizeof(tle_array_t));

  if (tle_array == NULL) {
    return NULL;
  }

  *tle_array = load_tles(filename);

  return tle_array;
}

void tle_free(void *tle_array) {
  free_tles(tle_array);
  free(tle_array);
}

// Setup and teardown functions for tests
int setup_nonexistent(void **state) {
  tle_array_t * tle_array = tle_load("tests/data/nonexistent.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int setup_empty(void **state) {
  tle_array_t * tle_array = tle_load("tests/data/empty.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int setup(void **state) {
  tle_array_t * tle_array = tle_load("tests/data/catalog.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int teardown(void **state) {
  tle_free(*state);

  return 0;
}

// Tests
void TLE_load_nonexistent_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 0);
}

void TLE_load_empty_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 0);
}

void TLE_load_invalid_index_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  tle_t * tle = get_orbit_by_index(&tle_array, 46);
  assert_null(tle);
}

void TLE_load_index_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  // AMS
  tle_t * tle = get_orbit_by_index(&tle_array, 9);
  assert_int_equal(tle->orbit.ep_year, 2023);
  assert_float_equal(tle->orbit.ep_day, 36.90397027, 1e-9);
  assert_float_equal(tle->orbit.rev, 15.27525761, 1e-9);
  assert_float_equal(tle->orbit.bstar, 0.0015803, 1e-9);
  assert_float_equal(tle->orbit.eqinc, 1.701910696, 1e-9);
  assert_float_equal(tle->orbit.ecc, 0.0022421, 1e-9);
  assert_float_equal(tle->orbit.mnan, 4.287516499, 1e-9);
  assert_float_equal(tle->orbit.argp, 2.001910105, 1e-9);
  assert_float_equal(tle->orbit.ascn, 2.717910487, 1e-9);
  // smjaxs is not set by sgdp4h
  assert_float_equal(tle->orbit.smjaxs, 0., 1e-9);
  assert_float_equal(tle->orbit.ndot2, 0.00042985, 1e-9);
  assert_float_equal(tle->orbit.nddot6, 0., 1e-9);
  assert_string_equal(tle->orbit.desig, "22057P");
  assert_int_equal(tle->orbit.norb, 3890);
  assert_int_equal(tle->orbit.satno, 52745);
}

void TLE_load_invalid_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  tle_t * tle = get_orbit_by_catalog_id(&tle_array, 12000);
  assert_null(tle);
}

void TLE_load_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  // ICEYE-X18
  tle_t * tle = get_orbit_by_catalog_id(&tle_array, 52749);
  assert_int_equal(tle->orbit.ep_year, 2023);
  assert_float_equal(tle->orbit.ep_day, 37.11421501, 1e-9);
  assert_float_equal(tle->orbit.rev, 15.15861349, 1e-9);
  assert_float_equal(tle->orbit.bstar, 0.00085173, 1e-9);
  assert_float_equal(tle->orbit.eqinc, 1.702179477, 1e-9);
  assert_float_equal(tle->orbit.ecc, 0.0009622, 1e-9);
  assert_float_equal(tle->orbit.mnan, 1.314424913, 1e-9);
  assert_float_equal(tle->orbit.argp, 0.888248671, 1e-9);
  assert_float_equal(tle->orbit.ascn, 2.685304246, 1e-9);
  // smjaxs is not set by sgdp4h
  assert_float_equal(tle->orbit.smjaxs, 0., 1e-9);
  assert_float_equal(tle->orbit.ndot2, 0.00016269, 1e-9);
  assert_float_equal(tle->orbit.nddot6, 0., 1e-9);
  assert_string_equal(tle->orbit.desig, "22057T");
  assert_int_equal(tle->orbit.norb, 3877);
  assert_int_equal(tle->orbit.satno, 52749);
}

// Entry point to run all tests
int run_tle_tests() {
  const struct CMUnitTest tests[] = {
    cmocka_unit_test_setup_teardown(TLE_load_nonexistent_file, setup_nonexistent, teardown),
    cmocka_unit_test_setup_teardown(TLE_load_empty_file, setup_empty, teardown),
    cmocka_unit_test_setup_teardown(TLE_load_invalid_index_from_file, setup, teardown),
    cmocka_unit_test_setup_teardown(TLE_load_index_from_file, setup, teardown),
    cmocka_unit_test_setup_teardown(TLE_load_invalid_catalog_id_from_file, setup, teardown),
    cmocka_unit_test_setup_teardown(TLE_load_catalog_id_from_file, setup, teardown),
  };

  return cmocka_run_group_tests_name("TLE", tests, NULL, NULL);
}
