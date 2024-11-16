#include "tests_rfio.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <cmocka.h>

#include "../rfio.h"

// Test sets
// There are 2 test sets consisting of 200 minutes of random 1024 S/s int16 IQ
// samples,
//
// they where generated as follow
// $ dd if=/dev/random of=random bs=12288000 count=4
// $ ./rffft -i random -p tests/bin/ -o t10_n120_32 -f 2250000000 -s 1024 -T 2016-09-28T00:20:00 -t 10 -n 120
// $ ./rffft -i random -p tests/bin/ -o t15_n90_8 -f 2250000000 -s 1024 -T 2016-09-28T00:20:00 -t 15 -n 90 -b
//
// Those 2 sets have different integration times, number of integrations per
// file and sample format The first one has 20 minutes per file, the second 22
// minutes and 30 seconds per file
//
// t10_n120_32 bin start times:
//  0: 00:20:00
//  1: 00:40:00
//  2: 01:00:00
//  3: 01:20:00
//  4: 01:40:00
//  5: 02:00:00
//  6: 02:20:00
//  7: 02:40:00
//  8: 03:00:00
//  9: 03:20:00
// 10: 03:40:00
//
// t15_n90_8 bin start times:
//  0: 00:20:00
//  1: 00:42:30
//  2: 01:05:00
//  3: 01:27:30
//  4: 01:50:00
//  5: 02:12:30
//  6: 02:35:00
//  7: 02:57:30
//  8: 03:20:00

// Tests

void IO_test_non_existent_file(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  start_bin = -1;
  num_integrations = -1;
  status =
      get_subs_from_datestrings("tests/bin/nonexistent", "2016-09-28T00:30:10",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, -1);
  assert_int_equal(start_bin, -1);
  assert_int_equal(num_integrations, -1);
}

void IO_test_only_start_time_int(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // Date before bin, so start at 0
  start_bin = -1;
  num_integrations = 3600;
  status =
      get_subs_from_datestrings("tests/bin/t10_n120_32", "2016-09-28T00:10:20",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 3600);

  // Start somewhere inside the bins aligned on start of bin
  start_bin = -1;
  num_integrations = 42;
  status =
      get_subs_from_datestrings("tests/bin/t10_n120_32", "2016-09-28T01:00:00",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 2);
  assert_int_equal(num_integrations, 42);

  // Start somewhere inside the bins
  start_bin = -1;
  num_integrations = 3;
  status =
      get_subs_from_datestrings("tests/bin/t10_n120_32", "2016-09-28T01:30:00",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 3);
  assert_int_equal(num_integrations, 3);

  start_bin = -1;
  num_integrations = 42;
  status =
      get_subs_from_datestrings("tests/bin/t10_n120_32", "2016-09-28T02:30:00",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 6);
  assert_int_equal(num_integrations, 42);
}

void IO_test_only_end_time_int(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // End somewhere inside the bins aligned on start of bin
  start_bin = 0;
  num_integrations = -1;
  status = get_subs_from_datestrings("tests/bin/t10_n120_32", NULL,
                                     "2016-09-28T01:00:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 360);

  // End somewhere inside the bins
  start_bin = 1;
  num_integrations = -1;
  status = get_subs_from_datestrings("tests/bin/t10_n120_32", NULL,
                                     "2016-09-28T01:04:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 1);
  assert_int_equal(num_integrations, 240);

  start_bin = 3;
  num_integrations = -1;
  status = get_subs_from_datestrings("tests/bin/t10_n120_32", NULL,
                                     "2016-09-28T01:30:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 3);
  assert_int_equal(num_integrations, 120);

  // End after last bin, update number of integrations to end of file
  start_bin = 6;
  num_integrations = 3600;
  status = get_subs_from_datestrings("tests/bin/t10_n120_32", NULL,
                                     "2016-09-28T04:30:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 6);
  assert_int_equal(num_integrations, 600);
}

void IO_test_only_start_time_char(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // Date before bin, so start at 0
  start_bin = -1;
  num_integrations = 3600;
  status =
      get_subs_from_datestrings("tests/bin/t15_n90_8", "2016-09-28T00:10:20",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 3600);

  // Start somewhere inside the bins aligned on start of bin
  start_bin = -1;
  num_integrations = 42;
  status =
      get_subs_from_datestrings("tests/bin/t15_n90_8", "2016-09-28T01:05:00",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 2);
  assert_int_equal(num_integrations, 42);

  // Start somewhere inside the bins
  start_bin = -1;
  num_integrations = 3;
  status =
      get_subs_from_datestrings("tests/bin/t15_n90_8", "2016-09-28T01:30:05",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 3);
  assert_int_equal(num_integrations, 3);

  start_bin = -1;
  num_integrations = 42;
  status =
      get_subs_from_datestrings("tests/bin/t15_n90_8", "2016-09-28T02:30:00",
                                NULL, &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 5);
  assert_int_equal(num_integrations, 42);
}

void IO_test_only_end_time_char(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // End somewhere inside the bins aligned on start of bin
  start_bin = 0;
  num_integrations = -1;
  status = get_subs_from_datestrings("tests/bin/t15_n90_8", NULL,
                                     "2016-09-28T01:05:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 270);

  // End somewhere inside the bins
  start_bin = 1;
  num_integrations = 42;
  status = get_subs_from_datestrings("tests/bin/t15_n90_8", NULL,
                                     "2016-09-28T01:04:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 1);
  assert_int_equal(num_integrations, 90);

  start_bin = 3;
  num_integrations = 0;
  status = get_subs_from_datestrings("tests/bin/t15_n90_8", NULL,
                                     "2016-09-28T01:30:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 3);
  assert_int_equal(num_integrations, 90);

  // End after last bin, update number of integrations to end of file
  start_bin = 6;
  num_integrations = 3600;
  status = get_subs_from_datestrings("tests/bin/t15_n90_8", NULL,
                                     "2016-09-28T04:30:00", &start_bin,
                                     &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 6);
  assert_int_equal(num_integrations, 270);
}

void IO_test_both_start_and_stop_int(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // Start before and end after bin files
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t10_n120_32", "2016-09-28T00:00:00", "2016-09-28T03:50:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 00);
  assert_int_equal(num_integrations, 1320);

  // Start and end corresponds to bin files
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t10_n120_32", "2016-09-28T00:20:00", "2016-09-28T03:40:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 1320);

  // Start and end inside bin
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t10_n120_32", "2016-09-28T00:30:00", "2016-09-28T03:30:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 1200);

  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t10_n120_32", "2016-09-28T00:40:00", "2016-09-28T03:20:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 1);
  assert_int_equal(num_integrations, 1080);

  // Start and end inverted inside bin
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t10_n120_32", "2016-09-28T03:30:00", "2016-09-28T00:30:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, -1);
  assert_int_equal(start_bin, 9);
  assert_int_equal(num_integrations, 100);
}

void IO_test_both_start_and_stop_char(void **state) {
  int status;
  int start_bin;
  int num_integrations;

  // Start before and end after bin files
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t15_n90_8", "2016-09-28T00:00:00", "2016-09-28T03:50:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 810);

  // Start and end corresponds to bin files
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t15_n90_8", "2016-09-28T00:20:00", "2016-09-28T03:20:00",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 810);

  // Start and end inside bin
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t15_n90_8", "2016-09-28T00:31:15", "2016-09-28T03:28:45",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 0);
  assert_int_equal(num_integrations, 810);

  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t15_n90_8", "2016-09-28T00:42:30", "2016-09-28T03:17:30",
      &start_bin, &num_integrations);
  assert_int_equal(status, 0);
  assert_int_equal(start_bin, 1);
  assert_int_equal(num_integrations, 630);

  // Start and end inverted inside bin
  start_bin = 10;
  num_integrations = 100;
  status = get_subs_from_datestrings(
      "tests/bin/t15_n90_8", "2016-09-28T03:28:45", "2016-09-28T00:31:15",
      &start_bin, &num_integrations);
  assert_int_equal(status, -1);
  assert_int_equal(start_bin, 8);
  assert_int_equal(num_integrations, 100);
}

// Entry point to run all tests
int run_io_tests() {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(IO_test_non_existent_file),
      cmocka_unit_test(IO_test_only_start_time_int),
      cmocka_unit_test(IO_test_only_end_time_int),
      cmocka_unit_test(IO_test_only_start_time_char),
      cmocka_unit_test(IO_test_only_end_time_char),
      cmocka_unit_test(IO_test_both_start_and_stop_int),
      cmocka_unit_test(IO_test_both_start_and_stop_char)
  };

  return cmocka_run_group_tests_name("IO", tests, NULL, NULL);
}
