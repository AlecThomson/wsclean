#include "fftkernels.h"

#include <cstring>

void fft2f_r2c_composite(size_t imgHeight, size_t imgWidth, const float *in,
                         fftwf_complex *out) {
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = imgHeight * complexWidth;

  fftwf_complex *temp1 = fftwf_alloc_complex(complexSize);

  fftwf_plan plan_r2c =
      fftwf_plan_dft_r2c_1d(imgWidth, nullptr, nullptr, FFTW_ESTIMATE);
  fftwf_plan plan_c2c = fftwf_plan_dft_1d(imgHeight, nullptr, nullptr,
                                          FFTW_FORWARD, FFTW_ESTIMATE);

  for (size_t y = 0; y < imgHeight; y++) {
    float *temp2 = fftwf_alloc_real(imgWidth);
    memcpy(temp2, &in[y * imgWidth], imgWidth * sizeof(float));
    fftwf_execute_dft_r2c(plan_r2c, temp2, &temp1[y * complexWidth]);
    fftwf_free(temp2);
  }

  for (size_t x = 0; x < complexWidth; x++) {
    fftwf_complex *temp2 = fftwf_alloc_complex(imgHeight);

    // Copy input
    for (size_t y = 0; y < imgHeight; y++) {
      memcpy(&temp2[y], &temp1[y * complexWidth + x], sizeof(fftwf_complex));
    }

    // Perform 1D FFT over column
    fftwf_execute_dft(plan_c2c, temp2, temp2);

    // Transpose output
    for (size_t y = 0; y < imgHeight; y++) {
      memcpy(&out[y * complexWidth + x], temp2[y], sizeof(fftwf_complex));
    }

    fftwf_free(temp2);
  }

  fftwf_free(temp1);

  fftwf_destroy_plan(plan_r2c);
  fftwf_destroy_plan(plan_c2c);
}

void fft2f_c2r_composite(size_t imgHeight, size_t imgWidth,
                         const fftwf_complex *in, float *out) {
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = imgHeight * complexWidth;

  fftwf_complex *temp1 = fftwf_alloc_complex(complexSize);

  fftwf_plan plan_c2c = fftwf_plan_dft_1d(imgHeight, nullptr, nullptr,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_plan plan_c2r =
      fftwf_plan_dft_c2r_1d(imgWidth, nullptr, nullptr, FFTW_ESTIMATE);

  for (size_t x = 0; x < complexWidth; x++) {
    // Transpose input
    for (size_t y = 0; y < imgHeight; y++) {
      memcpy(&temp1[x * imgHeight + y], &in[y * complexWidth + x],
             sizeof(fftwf_complex));
    }

    // Perform 1D C2C FFT over column
    fftwf_complex *temp1_ptr = &temp1[x * imgHeight];
    fftwf_execute_dft(plan_c2c, temp1_ptr, temp1_ptr);
  }

  for (size_t y = 0; y < imgHeight; y++) {
    fftwf_complex *temp2 = fftwf_alloc_complex(complexWidth);

    // Transpose input
    for (size_t x = 0; x < complexWidth; x++) {
      memcpy(&temp2[x], &temp1[x * imgHeight + y], sizeof(fftwf_complex));
    }

    // Perform 1D C2R FFT over row
    fftwf_execute_dft_c2r(plan_c2r, temp2, reinterpret_cast<float *>(temp2));

    // Copy output
    memcpy(&out[y * imgWidth], temp2, imgWidth * sizeof(float));

    fftwf_free(temp2);
  }

  fftwf_free(temp1);

  fftwf_destroy_plan(plan_c2c);
  fftwf_destroy_plan(plan_c2r);
}