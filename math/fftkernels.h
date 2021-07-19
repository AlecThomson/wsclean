#ifndef FFT_KERNELS_H
#define FFT_KERNELS_H

#include <fftw3.h>

void fft2f_r2c_composite(size_t imgHeight, size_t imgWidth, const float *in,
                         fftwf_complex *out, size_t threadCount = 1);

void fft2f_c2r_composite(size_t imgHeight, size_t imgWidth,
                         const fftwf_complex *in, float *out,
                         size_t threadCount = 1);

#endif