#include "multiscaletransforms.h"

#include <schaapcommon/fft/convolver.h>

using aocommon::Image;

void MultiScaleTransforms::Transform(std::vector<Image>& images, Image& scratch,
                                     float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  scratch = 0.0;

  schaapcommon::fft::Convolver::PrepareSmallKernel(
      scratch.Data(), _width, _height, shape.Data(), kernelSize, _threadCount);
  for (Image& image : images)
    schaapcommon::fft::Convolver::ConvolveSameSize(_fftwManager, image.Data(),
                                                   scratch.Data(), _width,
                                                   _height, _threadCount);
}

void MultiScaleTransforms::PrepareTransform(float* kernel, float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  std::fill_n(kernel, _width * _height, 0.0);

  schaapcommon::fft::Convolver::PrepareSmallKernel(
      kernel, _width, _height, shape.Data(), kernelSize, _threadCount);
}

void MultiScaleTransforms::FinishTransform(float* image, const float* kernel) {
  schaapcommon::fft::Convolver::ConvolveSameSize(_fftwManager, image, kernel,
                                                 _width, _height, _threadCount);
}
