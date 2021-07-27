#include "fftkernels.h"

#include <cstring>
#include <algorithm>

void fft2f_r2c_composite(fftwf_plan plan_r2c, fftwf_plan plan_c2c,
                         size_t imgHeight, size_t imgWidth, const float *in,
                         fftwf_complex *out,
                         aocommon::StaticFor<size_t> &loop) {
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = imgHeight * complexWidth;

  fftwf_complex *temp1 = fftwf_alloc_complex(complexSize);

  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    fftwf_complex *temp2 = fftwf_alloc_complex(complexWidth);
    float *temp2_ptr = reinterpret_cast<float *>(temp2);
    for (size_t y = yStart; y < yEnd; y++) {
      float *temp1_ptr = reinterpret_cast<float *>(&temp1[y * complexWidth]);
      std::copy_n(&in[y * imgWidth], imgWidth, temp2_ptr);
      fftwf_execute_dft_r2c(plan_r2c, temp2_ptr, temp2);
      std::copy_n(temp2_ptr, 2 * complexWidth, temp1_ptr);
    }
    fftwf_free(temp2);
  });

  loop.Run(0, complexWidth, [&](size_t xStart, size_t xEnd) {
    // Partially unroll over columns
    size_t unroll = 4;
    fftwf_complex *temp2 = fftwf_alloc_complex(unroll * imgHeight);

    for (size_t x = xStart; x < xEnd; x += unroll) {
      // Copy input
      for (size_t y = 0; y < imgHeight; y++) {
        for (size_t i = 0; i < unroll; i++) {
          if ((x + i) < xEnd) {
            float *temp1_ptr =
                reinterpret_cast<float *>(&temp1[y * complexWidth + x + i]);
            float *temp2_ptr =
                reinterpret_cast<float *>(&temp2[i * imgHeight + y]);
            std::copy_n(temp1_ptr, 2, temp2_ptr);
          }
        }
      }

      // Perform 1D FFT over columns
      for (size_t i = 0; i < unroll; i++) {
        fftwf_complex *temp2_ptr = &temp2[i * imgHeight];
        fftwf_execute_dft(plan_c2c, temp2_ptr, temp2_ptr);
      }

      // Transpose output
      for (size_t y = 0; y < imgHeight; y++) {
        for (size_t i = 0; i < unroll; i++) {
          if ((x + i) < xEnd) {
            float *temp2_ptr =
                reinterpret_cast<float *>(&temp2[i * imgHeight + y]);
            float *out_ptr =
                reinterpret_cast<float *>(&out[y * complexWidth + x + i]);
            std::copy_n(temp2_ptr, 2, out_ptr);
          }
        }
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
    size_t unroll = 4;
    for (size_t x = xStart; x < xEnd; x += unroll) {
      // Transpose input
      for (size_t y = 0; y < imgHeight; y++) {
        for (size_t i = 0; i < unroll; i++) {
          if ((x + i) < xEnd) {
            const float *in_ptr =
                reinterpret_cast<const float *>(&in[y * complexWidth + x + i]);
            float *temp1_ptr =
                reinterpret_cast<float *>(&temp1[(x + i) * imgHeight + y]);
            std::copy_n(in_ptr, 2, temp1_ptr);
          }
        }
      }

      // Perform 1D C2C FFT over columns
      for (size_t i = 0; i < unroll; i++) {
        if ((x + i) < xEnd) {
          fftwf_complex *temp1_ptr = &temp1[(x + i) * imgHeight];
          fftwf_execute_dft(plan_c2c, temp1_ptr, temp1_ptr);
        }
      }
    }
  });

  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    fftwf_complex *temp2 = fftwf_alloc_complex(complexWidth);
    for (size_t y = yStart; y < yEnd; y++) {
      // Transpose input
      for (size_t x = 0; x < complexWidth; x++) {
        float *temp1_ptr = reinterpret_cast<float *>(&temp1[x * imgHeight + y]);
        float *temp2_ptr = reinterpret_cast<float *>(&temp2[x]);
        std::copy_n(temp1_ptr, 2, temp2_ptr);
      }

      // Perform 1D C2R FFT over row
      fftwf_execute_dft_c2r(plan_c2r, temp2, reinterpret_cast<float *>(temp2));

      // Copy output
      float *temp2_ptr = reinterpret_cast<float *>(temp2);
      std::copy_n(temp2_ptr, imgWidth, &out[y * imgWidth]);
    }

    fftwf_free(temp2);
  });

  fftwf_free(temp1);
}