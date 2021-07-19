#include "fftconvolver.h"
#include "fftkernels.h"

#include "../system/fftwmanager.h"

#include <aocommon/uvector.h>
#include <aocommon/staticfor.h>

#include <fftw3.h>

#include <complex>
#include <stdexcept>

void FFTConvolver::Convolve(FFTWManager& fftw, float* image, size_t imgWidth,
                            size_t imgHeight, const float* kernel,
                            size_t kernelSize, size_t threadCount) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);
  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight,
                   threadCount);
}

void FFTConvolver::ReverseAndConvolve(class FFTWManager& fftw, float* image,
                                      size_t imgWidth, size_t imgHeight,
                                      const float* kernel, size_t kernelSize,
                                      size_t threadCount) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);

  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight,
                   threadCount);
}

void FFTConvolver::PrepareSmallKernel(float* dest, size_t imgWidth,
                                      size_t imgHeight, const float* kernel,
                                      size_t kernelSize) {
  if (kernelSize > imgWidth || kernelSize > imgHeight)
    throw std::runtime_error("Kernel size > image dimension");
  const float* kernelIter = kernel;
  for (size_t y = 0; y != kernelSize / 2; ++y) {
    size_t destY = imgHeight - kernelSize / 2 + y;
    size_t firstX = imgWidth - kernelSize / 2;
    float* destIter = &dest[destY * imgWidth + firstX];
    for (size_t x = 0; x != kernelSize / 2; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
    destIter = &dest[destY * imgWidth];
    for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
  }
  for (size_t y = kernelSize / 2; y != kernelSize; ++y) {
    size_t firstX = imgWidth - kernelSize / 2;
    float* destIter = &dest[firstX + (y - kernelSize / 2) * imgWidth];
    for (size_t x = 0; x != kernelSize / 2; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
    destIter = &dest[(y - kernelSize / 2) * imgWidth];
    for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
      *destIter = *kernelIter;
      ++kernelIter;
      ++destIter;
    }
  }
}

void FFTConvolver::PrepareKernel(float* dest, const float* source,
                                 size_t imgWidth, size_t imgHeight) {
  const float* sourceIter = source;
  for (size_t y = 0; y != imgHeight / 2; ++y) {
    size_t destY = imgHeight - imgHeight / 2 + y;
    size_t firstX = imgWidth - imgWidth / 2;
    float* destIter = &dest[destY * imgWidth + firstX];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
    destIter = &dest[destY * imgWidth];
    for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
  }
  for (size_t y = imgHeight / 2; y != imgHeight; ++y) {
    size_t firstX = imgWidth - imgWidth / 2;
    float* destIter = &dest[firstX + (y - imgHeight / 2) * imgWidth];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
    destIter = &dest[(y - imgHeight / 2) * imgWidth];
    for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
      *destIter = *sourceIter;
      ++sourceIter;
      ++destIter;
    }
  }
}

void FFTConvolver::ConvolveSameSize(FFTWManager& fftw, float* image,
                                    const float* kernel, size_t imgWidth,
                                    size_t imgHeight, size_t threadCount) {
  const size_t imgSize = imgWidth * imgHeight;
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = complexWidth * imgHeight;
  float* tempData = fftwf_alloc_real(imgSize);
  fftwf_complex* fftImageData = fftwf_alloc_complex(complexSize);
  fftwf_complex* fftKernelData = fftwf_alloc_complex(complexSize);

  std::unique_lock<std::mutex> lock(fftw.Mutex());
  fft2f_r2c_composite(imgHeight, imgWidth, image, fftImageData, threadCount);
  lock.unlock();

  std::copy_n(kernel, imgSize, tempData);
  lock.lock();
  fft2f_r2c_composite(imgHeight, imgWidth, tempData, fftKernelData);
  lock.unlock();

  float fact = 1.0 / imgSize;
  aocommon::StaticFor<size_t> loop(threadCount);
  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    for (size_t y = yStart; y != yEnd; ++y) {
      for (size_t x = 0; x != complexWidth; ++x) {
        size_t i = y * complexWidth + x;
        reinterpret_cast<std::complex<float>*>(fftImageData)[i] *=
            fact * reinterpret_cast<std::complex<float>*>(fftKernelData)[i];
      }
    }
  });

  lock.lock();
  fft2f_c2r_composite(imgHeight, imgWidth, fftImageData, image, threadCount);
  lock.unlock();

  fftwf_free(fftImageData);
  fftwf_free(fftKernelData);
  fftwf_free(tempData);
}

void FFTConvolver::Reverse(float* image, size_t imgWidth, size_t imgHeight) {
  for (size_t y = 0; y != imgHeight / 2; ++y) {
    size_t destY = imgHeight - 1 - y;
    float* sourcePtr = &image[y * imgWidth];
    float* destPtr = &image[destY * imgWidth];
    for (size_t x = 0; x != imgWidth / 2; ++x) {
      size_t destX = imgWidth - 1 - x;
      std::swap(sourcePtr[x], destPtr[destX]);
    }
  }
}
