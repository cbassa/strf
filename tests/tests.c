#include "tests_rffft_internal.h"
#include "tests_rftles.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmocka.h>

int main(void) {
  int failures = 0;

  failures += run_rffft_internal_tests();
  failures += run_tle_tests();

  return failures;
}
