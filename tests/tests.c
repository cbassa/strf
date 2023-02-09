#include "tests_rftles.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmocka.h>

int main(void) {
  int failures = 0;

  failures += run_tle_tests();

  return failures;
}
