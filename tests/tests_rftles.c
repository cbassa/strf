#include "tests_rftles.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <math.h>
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

int setup_alpha5(void **state) {
  tle_array_t * tle_array = load_tles("tests/data/alpha5.tle");

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
  assert_int_equal(tle_array.number_of_elements, 49);

  tle_t * tle = get_tle_by_index(&tle_array, 50);
  assert_null(tle);
}

void TLE_load_index_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 49);

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

  // 45th element, 60955 - 1998-067WV
  tle = get_tle_by_index(&tle_array, 45);
  assert_int_equal(tle->orbit.satno, 60955);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "1998-067WV");

  // 46th element, 45727 - 2 PATHFINDER
  tle = get_tle_by_index(&tle_array, 46);
  assert_int_equal(tle->orbit.satno, 45727);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "2 PATHFINDER");

  // 47th element, 5485 - SHINSEI (MS-F2)
  tle = get_tle_by_index(&tle_array, 47);
  assert_int_equal(tle->orbit.satno, 5485);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "SHINSEI (MS-F2)");

  // 48th element, 07530 - OSCAR 7 (AO-7)
  tle = get_tle_by_index(&tle_array, 48);
  assert_int_equal(tle->orbit.satno, 7530);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "OSCAR 7 (AO-7)");
}

void TLE_load_invalid_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 49);

  tle_t * tle = get_tle_by_catalog_id(&tle_array, 12000);
  assert_null(tle);
}

void TLE_load_catalog_id_from_file(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 49);

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

  // 45th element, 60955 - 1998-067WV
  tle = get_tle_by_catalog_id(&tle_array, 60955);
  assert_int_equal(tle->orbit.satno, 60955);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "1998-067WV");

  // 46th element, 45727 - 2 PATHFINDER
  tle = get_tle_by_catalog_id(&tle_array, 45727);
  assert_int_equal(tle->orbit.satno, 45727);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "2 PATHFINDER");

  // 47th element, 5485 - SHINSEI (MS-F2)
  tle = get_tle_by_catalog_id(&tle_array, 5485);
  assert_int_equal(tle->orbit.satno, 5485);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "SHINSEI (MS-F2)");

  // 48th element, 07530 - OSCAR 7 (AO-7)
  tle = get_tle_by_catalog_id(&tle_array, 7530);
  assert_int_equal(tle->orbit.satno, 7530);
  assert_non_null(tle->name);
  assert_string_equal(tle->name, "OSCAR 7 (AO-7)");
}

void TLE_decode_alpha5_designation(void **state) {
  tle_array_t tle_array = **(tle_array_t **)state;

  assert_non_null(tle_array.tles);
  assert_int_equal(tle_array.number_of_elements, 2);

  // 1st element, B5544 - ISS
  tle_t * tle = get_tle_by_catalog_id(&tle_array, 115544);

  assert_non_null(tle->name);
  assert_string_equal(tle->name, "ISS (ZARYA)");

  assert_string_equal(tle->orbit.desig, "98067A");
  assert_int_equal(tle->orbit.satno, 115544);

  // 2nd element, 339999 - ISS
  tle = get_tle_by_catalog_id(&tle_array, 339999);

  assert_non_null(tle->name);
  assert_string_equal(tle->name, "ISS (ZARYA)");

  assert_string_equal(tle->orbit.desig, "98067A");
  assert_int_equal(tle->orbit.satno, 339999);

}

// The same four objects in every element set format (tests/data/omm.*):
// ISS, GOES 18, MERIDIAN 7 and CREW DRAGON 13 (100882, alpha-5 A0882 in
// the TLE file). The OMM encodings must give identical elements, the TLE
// file the same elements within TLE precision.
static const char *omm_files[] = {
  "tests/data/omm.json", "tests/data/omm.xml", "tests/data/omm.kvn", "tests/data/omm.csv"
};
static const long omm_satnos[] = {25544, 51850, 40296, 100882};
static const char *omm_names[] = {"ISS (ZARYA)", "GOES 18", "MERIDIAN 7", "CREW DRAGON 13"};

static void assert_same_orbit(const orbit_t *a, const orbit_t *b) {
  assert_int_equal(a->satno, b->satno);
  assert_int_equal(a->ep_year, b->ep_year);
  assert_true(a->ep_day == b->ep_day);
  assert_true(a->rev == b->rev);
  assert_true(a->bstar == b->bstar);
  assert_true(a->eqinc == b->eqinc);
  assert_true(a->ecc == b->ecc);
  assert_true(a->mnan == b->mnan);
  assert_true(a->argp == b->argp);
  assert_true(a->ascn == b->ascn);
  assert_true(a->ndot2 == b->ndot2);
  assert_true(a->nddot6 == b->nddot6);
  assert_int_equal(a->norb, b->norb);
  assert_string_equal(a->desig, b->desig);
}

void TLE_load_omm_encodings(void **state) {
  (void)state;
  tle_array_t *ref = load_tles((char *)omm_files[0]);

  assert_non_null(ref);
  assert_int_equal(ref->number_of_elements, 4);
  for (int i = 0; i < 4; i++) {
    assert_int_equal(ref->tles[i].orbit.satno, omm_satnos[i]);
    assert_non_null(ref->tles[i].name);
    assert_string_equal(ref->tles[i].name, omm_names[i]);
  }

  // ISS, field by field
  tle_t *tle = get_tle_by_catalog_id(ref, 25544);
  assert_non_null(tle);
  assert_int_equal(tle->orbit.ep_year, 2026);
  assert_float_equal(tle->orbit.ep_day, 277.07866052, 1e-8);
  assert_float_equal(tle->orbit.rev, 15.48727366, 1e-12);
  assert_float_equal(tle->orbit.bstar, 0.71811468e-4, 1e-15);
  assert_float_equal(tle->orbit.eqinc, 51.6316 * M_PI / 180.0, 1e-12);
  assert_float_equal(tle->orbit.ecc, 0.0006887, 1e-12);
  assert_float_equal(tle->orbit.ndot2, 0.00003472, 1e-15);
  assert_string_equal(tle->orbit.desig, "98067A");
  assert_int_equal(tle->orbit.norb, 58863);

  for (int f = 1; f < 4; f++) {
    tle_array_t *other = load_tles((char *)omm_files[f]);
    assert_non_null(other);
    assert_int_equal(other->number_of_elements, 4);
    for (int i = 0; i < 4; i++) {
      assert_same_orbit(&other->tles[i].orbit, &ref->tles[i].orbit);
      assert_string_equal(other->tles[i].name, ref->tles[i].name);
    }
    free_tles(other);
  }

  free_tles(ref);
}

void TLE_load_tle_and_omm(void **state) {
  (void)state;
  tle_array_t *omm = load_tles("tests/data/omm.json");
  tle_array_t *tle = load_tles("tests/data/omm.tle");
  const double ang = 1e-4 * M_PI / 180.0 * 1.001;

  assert_non_null(omm);
  assert_non_null(tle);
  assert_int_equal(tle->number_of_elements, 4);

  for (int i = 0; i < 4; i++) {
    const orbit_t *o = &omm->tles[i].orbit, *t = &tle->tles[i].orbit;

    assert_int_equal(t->satno, omm_satnos[i]);
    assert_string_equal(tle->tles[i].name, omm->tles[i].name);
    assert_string_equal(t->desig, o->desig);
    assert_int_equal(t->ep_year, o->ep_year);
    assert_float_equal(t->ep_day, o->ep_day, 1e-8);
    assert_float_equal(t->rev, o->rev, 1e-8);
    assert_float_equal(t->ecc, o->ecc, 1e-7);
    assert_float_equal(t->eqinc, o->eqinc, ang);
    assert_float_equal(t->ascn, o->ascn, ang);
    assert_float_equal(t->argp, o->argp, ang);
    assert_float_equal(t->mnan, o->mnan, ang);
    assert_float_equal(t->bstar, o->bstar, 1e-4 * fabs(o->bstar) + 1e-14);
    assert_float_equal(t->ndot2, o->ndot2, 1e-8);
  }

  free_tles(omm);
  free_tles(tle);
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
    cmocka_unit_test_setup_teardown(TLE_decode_alpha5_designation, setup_alpha5, teardown),    
    cmocka_unit_test(TLE_load_omm_encodings),
    cmocka_unit_test(TLE_load_tle_and_omm),
  };

  return cmocka_run_group_tests_name("TLE", tests, NULL, NULL);
}
