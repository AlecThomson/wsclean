#include "multiscaletransforms.h"

#include <schaapcommon/fft/fftconvolver.h>

using aocommon::Image;
using schaapcommon::fft::FftConvolver;

void MultiScaleTransforms::Transform(std::vector<Image>& images, Image& scratch,
                                     float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  scratch = 0.0;

  FftConvolver::PrepareSmallKernel(scratch.Data(), _width, _height,
                                   shape.Data(), kernelSize, _threadCount);
  for (Image& image : images)
    FftConvolver::ConvolveSameSize(_fftwManager, image.Data(), scratch.Data(),
                                   _width, _height, _threadCount);
}

void MultiScaleTransforms::PrepareTransform(float* kernel, float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  std::fill_n(kernel, _width * _height, 0.0);

  FftConvolver::PrepareSmallKernel(kernel, _width, _height, shape.Data(),
                                   kernelSize, _threadCount);
}

void MultiScaleTransforms::FinishTransform(float* image, const float* kernel) {
  FftConvolver::ConvolveSameSize(_fftwManager, image, kernel, _width, _height,
                                 _threadCount);
}
