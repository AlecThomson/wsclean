#include "fftkernels.h"

#include <cstring>

void fft2f_r2c_composite(fftwf_plan plan_r2c, fftwf_plan plan_c2c,
                         size_t imgHeight, size_t imgWidth, const float *in,
                         fftwf_complex *out,
                         aocommon::StaticFor<size_t> &loop) {
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = imgHeight * complexWidth;

  fftwf_complex *temp1 = fftwf_alloc_complex(complexSize);

  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    fftwf_complex *temp2 = fftwf_alloc_complex(complexWidth);
    for (size_t y = yStart; y < yEnd; y++) {
      memcpy(temp2, &in[y * imgWidth], imgWidth * sizeof(float));
      fftwf_execute_dft_r2c(plan_r2c, reinterpret_cast<float *>(temp2), temp2);
      memcpy(&temp1[y * complexWidth], temp2,
             complexWidth * sizeof(fftwf_complex));
    }
    fftwf_free(temp2);
  });

  loop.Run(0, complexWidth, [&](size_t xStart, size_t xEnd) {
    fftwf_complex *temp2 = fftwf_alloc_complex(imgHeight);

    for (size_t x = xStart; x < xEnd; x++) {
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
    }

    fftwf_free(temp2);
  });

  fftwf_free(temp1);
}

void fft2f_c2r_composite(fftwf_plan plan_c2c, fftwf_plan plan_c2r,
                         size_t imgHeight, size_t imgWidth,
                         const fftwf_complex *in, float *out,
                         aocommon::StaticFor<size_t> &loop) {
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = imgHeight * complexWidth;

  fftwf_complex *temp1 = fftwf_alloc_complex(complexSize);

  loop.Run(0, complexWidth, [&](size_t xStart, size_t xEnd) {
    for (size_t x = xStart; x < xEnd; x++) {
      // Transpose input
      for (size_t y = 0; y < imgHeight; y++) {
        memcpy(&temp1[x * imgHeight + y], &in[y * complexWidth + x],
               sizeof(fftwf_complex));
      }

      // Perform 1D C2C FFT over column
      fftwf_complex *temp1_ptr = &temp1[x * imgHeight];
      fftwf_execute_dft(plan_c2c, temp1_ptr, temp1_ptr);
    }
  });

  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    fftwf_complex *temp2 = fftwf_alloc_complex(complexWidth);
    for (size_t y = yStart; y < yEnd; y++) {
      // Transpose input
      for (size_t x = 0; x < complexWidth; x++) {
        memcpy(&temp2[x], &temp1[x * imgHeight + y], sizeof(fftwf_complex));
      }

      // Perform 1D C2R FFT over row
      fftwf_execute_dft_c2r(plan_c2r, temp2, reinterpret_cast<float *>(temp2));

      // Copy output
      memcpy(&out[y * imgWidth], temp2, imgWidth * sizeof(float));
    }

    fftwf_free(temp2);
  });

  fftwf_free(temp1);
}