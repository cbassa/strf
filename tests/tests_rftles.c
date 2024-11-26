#include "tests_rftles.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmocka.h>

#include "../rftles.h"

// Setup and teardown functions for tests
int setup_nonexistent(void **state) {
  tle_array_t * tle_array = load_tles("tests/data/nonexistent.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int setup_empty(void **state) {
  tle_array_t * tle_array = load_tles("tests/data/empty.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int setup(void **state) {
  tle_array_t * tle_array = load_tles("tests/data/catalog.tle");

  if (tle_array == NULL) {
    return -1;
  }

  *state = tle_array;

  return 0;
}

int teardown(void **state) {
  free_tles(*state);

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

  tle_t * tle = get_tle_by_index(&tle_array, 46);
  assert_null(tle);
}

void TLE_load_index_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  // 9th element, 52745 - AMS
  tle_t * tle = get_tle_by_index(&tle_array, 9);

  assert_non_null(tle->name);
  assert_string_equal(tle->name, "AMS");

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

  // 7th element, 52743 - OBJECT M, no name
  tle = get_tle_by_index(&tle_array, 7);
  assert_null(tle->name);

  // 0th element, 52736 - LEMUR 2 KAREN_B
  tle = get_tle_by_index(&tle_array, 0);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "LEMUR 2 KAREN_B");
}

void TLE_load_invalid_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  tle_t * tle = get_tle_by_catalog_id(&tle_array, 12000);
  assert_null(tle);
}

void TLE_load_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 45);

  // 13th element, 52749, - ICEYE-X18
  tle_t * tle = get_tle_by_catalog_id(&tle_array, 52749);

  assert_non_null(tle->name);
  assert_string_equal(tle->name, "ICEYE-X18");

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

  // 7th element, 52743 - OBJECT M, no name
  tle = get_tle_by_catalog_id(&tle_array, 52743);
  assert_null(tle->name);

  // 0th element, 52736 - LEMUR 2 KAREN_B
  tle = get_tle_by_catalog_id(&tle_array, 52736);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "LEMUR 2 KAREN_B");
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
