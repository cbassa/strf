#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "rftime.h"
#include "rfio.h"

#define MAX_REJECT 0.5
#define MIN_NPIXELS 5
#define GOOD_PIXEL 0
#define BAD_PIXEL 1
#define KREJ 2.5
#define MAX_ITERATIONS 5

void zsc_sample(struct spectrogram *image, int maxpix, double *samples, int *nsamples) {
    int nc = image->nchan;
    int nl = image->nsub;
    double stride = fmax(1.0, sqrt((nc - 1) * (nl - 1) / (double)maxpix));
    int istride = (int)stride;

    int count = 0;
    for (int i = 0; i < nc; i += istride) {
        for (int j = 0; j < nl; j += istride) {
            if (count >= maxpix) {
                *nsamples = count;
                return;
            }
            samples[count++] = image->z[i * nl + j];
        }
    }
    *nsamples = count;
}

void zsc_compute_sigma(double *flat, int *badpix, int npix, int *ngoodpix, double *mean, double *sigma) {
    double sumz = 0.0;
    double sumsq = 0.0;
    *ngoodpix = 0;

    for (int i = 0; i < npix; i++) {
        if (badpix[i] == GOOD_PIXEL) {
            sumz += flat[i];
            sumsq += flat[i] * flat[i];
            (*ngoodpix)++;
        }
    }

    if (*ngoodpix == 0) {
        *mean = NAN;
        *sigma = NAN;
    } else if (*ngoodpix == 1) {
        *mean = sumz;
        *sigma = NAN;
    } else {
        *mean = sumz / *ngoodpix;
        double temp = sumsq / (*ngoodpix - 1) - (sumz * sumz) / (*ngoodpix * (*ngoodpix - 1));
        *sigma = temp < 0.0 ? 0.0 : sqrt(temp);
    }
}

void zsc_fit_line(double *samples, int npix, double krej, int ngrow, int maxiter,
                  int *ngoodpix_out, double *zstart, double *zslope) {
    double xscale = 2.0 / (npix - 1);
    double *xnorm = (double *)malloc(npix * sizeof(double));
    int *badpix = (int *)calloc(npix, sizeof(int));
    int *conv = (int *)calloc(npix, sizeof(int));

    for (int i = 0; i < npix; i++) {
        xnorm[i] = i * xscale - 1.0;
    }

    int ngoodpix = npix;
    int minpix = fmax(MIN_NPIXELS, (int)(npix * MAX_REJECT));
    int last_ngoodpix = npix + 1;

    double intercept = 0.0;
    double slope = 0.0;

    for (int niter = 0; niter < maxiter; niter++) {
        if (ngoodpix >= last_ngoodpix || ngoodpix < minpix) {
            break;
        }

        double sumx = 0.0, sumxx = 0.0, sumxy = 0.0, sumy = 0.0;
        int count = 0;
        for (int i = 0; i < npix; i++) {
            if (badpix[i] == GOOD_PIXEL) {
                sumx += xnorm[i];
                sumxx += xnorm[i] * xnorm[i];
                sumxy += xnorm[i] * samples[i];
                sumy += samples[i];
                count++;
            }
        }

        double delta = count * sumxx - sumx * sumx;
        intercept = (sumxx * sumy - sumx * sumxy) / delta;
        slope = (count * sumxy - sumx * sumy) / delta;

        double *fitted = (double *)malloc(npix * sizeof(double));
        double *flat = (double *)malloc(npix * sizeof(double));
        for (int i = 0; i < npix; i++) {
            fitted[i] = xnorm[i] * slope + intercept;
            flat[i] = samples[i] - fitted[i];
        }

        double mean, sigma;
        zsc_compute_sigma(flat, badpix, npix, &ngoodpix, &mean, &sigma);

        double threshold = sigma * krej;
        double lcut = -threshold;
        double hcut = threshold;

        for (int i = 0; i < npix; i++) {
            if (flat[i] < lcut || flat[i] > hcut) {
                badpix[i] = BAD_PIXEL;
            }
        }
        // temp = np.convolve(badpix, np.ones(ngrow), mode='same')
        // fixed convolution computation
        for (int i = 0; i < npix; i++) {
            int sum = 0;
            for (int j = -ngrow/2; j < ngrow/2; j++) {
                int idx = i + j;
                if (idx >= 0 && idx < npix) {
                    sum+= badpix[idx];
                }
            }
            conv[i] = sum;
        }

        free(fitted);
        free(flat);

        last_ngoodpix = ngoodpix;
        ngoodpix = 0;
        for (int i = 0; i < npix; i++) {
            badpix[i] = conv[i];
            if (conv[i] == GOOD_PIXEL) {
                ngoodpix++;
            }
        }
    }

    *zstart = intercept - slope;
    *zslope = slope * xscale;
    *ngoodpix_out = ngoodpix; 
    free(xnorm);
    free(badpix);
    free(conv);
}

int compare_doubles(const void *a, const void *b) {
    double diff = (*(double *)a - *(double *)b);
    return (diff > 0) - (diff < 0);
}

void zscale(struct spectrogram *image, int nsamples, double contrast, double *z1, double *z2) {
    double *samples = (double *)malloc(nsamples * sizeof(double));
    int npix;

    zsc_sample(image, nsamples, samples, &npix);

    qsort(samples, npix, sizeof(double), compare_doubles);

    double zmin = samples[0];
    double zmax = samples[npix - 1];
    double median;
    int center_pixel = (npix - 1) / 2;
    if (npix % 2 == 1) {
        median = samples[center_pixel];
    } else {
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1]);
    }

    int ngoodpix;
    double zstart, zslope;
    zsc_fit_line(samples, npix, KREJ, fmax(1, npix * 0.01), MAX_ITERATIONS,
                 &ngoodpix, &zstart, &zslope);


    if (ngoodpix < fmax(MIN_NPIXELS, npix * MAX_REJECT)) {
        *z1 = zmin;
        *z2 = zmax;
    } else {
        if (contrast > 0) zslope /= contrast;
        *z1 = fmax(zmin, median - (center_pixel - 1) * zslope);
        *z2 = fmin(zmax, median + (npix - center_pixel) * zslope);
    }

    free(samples);
}
