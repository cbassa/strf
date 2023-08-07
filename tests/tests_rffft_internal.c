#include "tests_rffft_internal.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmocka.h>

#include "../rffft_internal.h"

// Tests

// Test SatDump filenames
void rffft_internal_parse_satdump_filenames(void **state) {
  double samplerate = 0;
  double frequency = 0;
  char format = '\0';
  char starttime[] = "YYYY-mm-ddTHH:MM:SS.sss";
  char ref_format = '\0';

  // s8 file without milliseconds
  ref_format = 'c';
  assert_int_equal(0, rffft_params_from_filename("2023-08-05_08-02-00_16000000SPS_2274000000Hz.s8", &samplerate, &frequency, &format, starttime));
  // assert_double_equal has been introduced in cmocka 1.1.6 not available on most distribs yet
  assert_float_equal(16e6, samplerate, 1e-12);
  assert_float_equal(2.274e9, frequency, 1e-12);
  assert_memory_equal(&ref_format, &format, 1);
  assert_string_equal("2023-08-05T08:02:00", starttime);

  // f32 file with milliseconds
  ref_format = 'f';
  assert_int_equal(0, rffft_params_from_filename("2023-08-17_11-41-14-373_1000000SPS_100000000Hz.f32", &samplerate, &frequency, &format, starttime));
  assert_float_equal(1e6, samplerate, 1e-12);
  assert_float_equal(100e6, frequency, 1e-12);
  assert_memory_equal(&ref_format, &format, 1);
  assert_string_equal("2023-08-17T11:41:14.373", starttime);

  // s16 file with broken milliseconds format
  ref_format = 'i';
  assert_int_equal(0, rffft_params_from_filename("2023-08-05_18-02-45-1691258565.534000_8000000SPS_2284000000Hz.s16", &samplerate, &frequency, &format, starttime));
  assert_float_equal(8e6, samplerate, 1e-12);
  assert_float_equal(2.284e9, frequency, 1e-12);
  assert_memory_equal(&ref_format, &format, 1);
  assert_string_equal("2023-08-05T18:02:45.534", starttime);

  // wav file with milliseconds
  ref_format = 'w';
  assert_int_equal(0, rffft_params_from_filename("2023-08-07_16-36-47-749_2400000SPS_100000000Hz.wav", &samplerate, &frequency, &format, starttime));
  assert_float_equal(2.4e6, samplerate, 1e-12);
  assert_float_equal(100e6, frequency, 1e-12);
  assert_memory_equal(&ref_format, &format, 1);
  assert_string_equal("2023-08-07T16:36:47.749", starttime);

  assert_int_equal(-1, rffft_params_from_filename("2023-08-05-19:59:30_16000000SPS_402000000Hz.f32", &samplerate, &frequency, &format, starttime));
}

// Test GQRX filenames
void rffft_internal_parse_gqrx_filenames(void **state) {
  double samplerate = 0;
  double frequency = 0;
  char format = '\0';
  char starttime[] = "YYYY-mm-ddTHH:MM:SS.sss";
  char ref_format = '\0';

  ref_format = 'f';
  assert_int_equal(0, rffft_params_from_filename("gqrx_20230806_151838_428000000_200000_fc.raw", &samplerate, &frequency, &format, starttime));
  // assert_double_equal has been introduced in cmocka 1.1.6 not available on most distribs yet
  assert_float_equal(200e3, samplerate, 1e-12);
  assert_float_equal(428e6, frequency, 1e-12);
  assert_memory_equal(&ref_format, &format, 1);
  assert_string_equal("2023-08-06T15:18:38", starttime);

  assert_int_equal(-1, rffft_params_from_filename("gqrx_2023-08-06_15:18:38_428000000_200000_fc.raw", &samplerate, &frequency, &format, starttime));
}

// Entry point to run all tests
int run_rffft_internal_tests() {
  const struct CMUnitTest tests[] = {
    cmocka_unit_test(rffft_internal_parse_satdump_filenames),
    cmocka_unit_test(rffft_internal_parse_gqrx_filenames),
  };

  return cmocka_run_group_tests_name("rffft internal", tests, NULL, NULL);
}
