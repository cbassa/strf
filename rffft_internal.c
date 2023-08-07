#include "rffft_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>


// Filename formats:
// - SatDump
//   - 2023-08-05_08-02-00_16000000SPS_2274000000Hz.s8
//   - 2023-08-05_18-02-45-534_16000000SPS_2284000000Hz.s16
//   - 2023-08-05_18-02-45-1691258565.534000_16000000SPS_2284000000Hz.f32
//   - 2023-08-07_16-36-47-1691426207.749000_2400000SPS_100000000Hz.wav
//   s8: char, s16 short int, f32 float.
//   SatDump also supports a compressed versions of s8/s16/f32 with .ziq
//   extension. Those are not yet supported
//   timestamp can have an added milliseconds field, configurable. This feature
//   was broken during some time so files with this convention still exists.
// - GQRX:
//   - gqrx_20230806_151838_428000000_200000_fc.raw
//   format always float32
int rffft_params_from_filename(char * filename, double * samplerate, double * frequency, char * format, char * starttime) {
  // Temp vars to hold parsed values
  int p_year, p_month, p_day, p_hours, p_minutes, p_seconds, p_fractal_seconds;
  int p_dummy_int;
  double p_samplerate, p_frequency;
  char p_format[16];
  char p_dummy_string[128];
  int parsed_tokens;

  char * base_filename = basename(filename);

  // Broken SatDump string with milliseconds activated
  parsed_tokens = sscanf(
    base_filename,
    "%04d-%02d-%02d_%02d-%02d-%02d-%d.%06d_%lfSPS_%lfHz.%s",
    &p_year,
    &p_month,
    &p_day,
    &p_hours,
    &p_minutes,
    &p_seconds,
    &p_dummy_int,
    &p_fractal_seconds,
    &p_samplerate,
    &p_frequency,
    p_format);

  if (parsed_tokens == 11) {
    *samplerate = p_samplerate;
    *frequency = p_frequency;

    if ((strlen(p_format) == 2) && (strncmp("s8", p_format, 2) == 0)) {
      *format = 'c';
    } else if ((strlen(p_format) == 3) && (strncmp("s16", p_format, 3) == 0)) {
      *format = 'i';
    } else if ((strlen(p_format) == 3) && (strncmp("f32", p_format, 3) == 0)) {
      *format = 'f';
    } else if ((strlen(p_format) == 3) && (strncmp("wav", p_format, 3) == 0)) {
      *format = 'w';
    } else {
      printf("Unsupported SatDump format %s\n", p_format);
      return -1;
    }

    snprintf(starttime, 32, "%04d-%02d-%02dT%02d:%02d:%02d.%03d", p_year, p_month, p_day, p_hours, p_minutes, p_seconds, p_fractal_seconds / 1000);

    return 0;
  }

  // SatDump string with milliseconds activated
  parsed_tokens = sscanf(
    base_filename,
    "%04d-%02d-%02d_%02d-%02d-%02d-%03d_%lfSPS_%lfHz.%s",
    &p_year,
    &p_month,
    &p_day,
    &p_hours,
    &p_minutes,
    &p_seconds,
    &p_fractal_seconds,
    &p_samplerate,
    &p_frequency,
    p_format);

  if (parsed_tokens == 10) {
    *samplerate = p_samplerate;
    *frequency = p_frequency;

    if ((strlen(p_format) == 2) && (strncmp("s8", p_format, 2) == 0)) {
      *format = 'c';
    } else if ((strlen(p_format) == 3) && (strncmp("s16", p_format, 3) == 0)) {
      *format = 'i';
    } else if ((strlen(p_format) == 3) && (strncmp("f32", p_format, 3) == 0)) {
      *format = 'f';
    } else if ((strlen(p_format) == 3) && (strncmp("wav", p_format, 3) == 0)) {
      *format = 'w';
    } else {
      printf("Unsupported SatDump format %s\n", p_format);
      return -1;
    }

    snprintf(starttime, 32, "%04d-%02d-%02dT%02d:%02d:%02d.%03d", p_year, p_month, p_day, p_hours, p_minutes, p_seconds, p_fractal_seconds);

    return 0;
  }

  // SatDump string without milliseconds activated
  parsed_tokens = sscanf(
    base_filename,
    "%04d-%02d-%02d_%02d-%02d-%02d_%lfSPS_%lfHz.%s",
    &p_year,
    &p_month,
    &p_day,
    &p_hours,
    &p_minutes,
    &p_seconds,
    &p_samplerate,
    &p_frequency,
    p_format);

  if (parsed_tokens == 9) {
    *samplerate = p_samplerate;
    *frequency = p_frequency;

    if ((strlen(p_format) == 2) && (strncmp("s8", p_format, 2) == 0)) {
      *format = 'c';
    } else if ((strlen(p_format) == 3) && (strncmp("s16", p_format, 3) == 0)) {
      *format = 'i';
    } else if ((strlen(p_format) == 3) && (strncmp("f32", p_format, 3) == 0)) {
      *format = 'f';
    } else if ((strlen(p_format) == 3) && (strncmp("wav", p_format, 3) == 0)) {
      *format = 'w';
    } else {
      printf("Unsupported SatDump format %s\n", p_format);
      return -1;
    }

    snprintf(starttime, 32, "%04d-%02d-%02dT%02d:%02d:%02d", p_year, p_month, p_day, p_hours, p_minutes, p_seconds);

    return 0;
  }

  // GQRX
  parsed_tokens = sscanf(
    base_filename,
    "gqrx_%04d%02d%02d_%02d%02d%02d_%lf_%lf_%s.raw",
    &p_year,
    &p_month,
    &p_day,
    &p_hours,
    &p_minutes,
    &p_seconds,
    &p_frequency,
    &p_samplerate,
    p_dummy_string);

  if (parsed_tokens == 9) {
    *samplerate = p_samplerate;
    *frequency = p_frequency;
    *format = 'f';

    snprintf(starttime, 32, "%04d-%02d-%02dT%02d:%02d:%02d", p_year, p_month, p_day, p_hours, p_minutes, p_seconds);

    return 0;
  }

  return -1;
}
