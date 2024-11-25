#ifndef _RFFFT_INTERNAL_H
#define _RFFFT_INTERNAL_H

#include "sgdp4h.h"

#ifdef __cplusplus
extern "C" {
#endif

// input:
// filename: filename string to parse
// output:
// samplerate: parsed samplerate
// frequency: parsed frequency
// format: parsed sample format: char: 'c', int: 'i', float: 'f', wav: 'w'
// starttime: parsed start time string formatted YYYY-MM-DDTHH:MM:SS.sss
int rffft_params_from_filename(char * filename, double * samplerate, double * frequency, char * format, char * starttime);

#ifdef __cplusplus
}
#endif

#endif /* _RFFFT_INTERNAL_H */
