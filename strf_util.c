/* > strf_util.c
 *
 * strf-specific pieces on top of the SGP4 library in extern/dundee_sgp4.
 */

#include "sgdp4h.h"

#include <stdarg.h>

/* The SGP4 library reports errors through fatal_error(), and its own
 * version exits. strf reports the error and carries on (for example when
 * one element set out of a catalogue fails to initialise), so it keeps its
 * own definition here; with this one linked, the library's ferror.o is not
 * pulled in. */
void fatal_error(const char *format, ...)
{
va_list arg_ptr;

    fflush(stdout);

    fprintf(stderr, "\nDundee Satellite Lab fatal run-time error:\n");

    va_start(arg_ptr, format);
    vfprintf(stderr, format, arg_ptr);
    va_end(arg_ptr);

    return;
}
