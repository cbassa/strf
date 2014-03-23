/* > satutl.c
 *
 */


#include "sgdp4h.h"

#include <stdarg.h>

void fatal_error(const char *format, ...)
{
va_list arg_ptr;

    fflush(stdout);

    fprintf(stderr, "\nDundee Satellite Lab fatal run-time error:\n");

    va_start(arg_ptr, format);
    vfprintf(stderr, format, arg_ptr);
    va_end(arg_ptr);

    fprintf(stderr, "\nNow terminating the program...\n");
    fflush(stderr);

    exit(5);

}

/* ===================================================================== */
