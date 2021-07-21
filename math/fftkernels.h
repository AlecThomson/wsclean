#ifndef FFT_KERNELS_H
#define FFT_KERNELS_H

#include <fftw3.h>

void fft2f_r2c_composite(fftwf_plan plan_r2c, fftwf_plan plan_c2c,
                         size_t imgHeight, size_t imgWidth, const float *in,
                         fftwf_complex *out);

void fft2f_c2r_composite(fftwf_plan plan_c2c, fftwf_plan plan_c2r,
                         size_t imgHeight, size_t imgWidth,
                         const fftwf_complex *in, float *out);

#endif