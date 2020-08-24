#include "wgriddinggridder_simple.h"
#include "gridder_cxx.h"

using namespace std;
using namespace ducc0;

WGriddingGridder_Simple::WGriddingGridder_Simple(
    size_t width, size_t height, size_t width_t, size_t height_t,
    double pixelSizeX, double pixelSizeY, size_t nthreads, double epsilon,
    size_t verbosity)
    : width_(width),
      height_(height),
      width_t_(width_t),
      height_t_(height_t),
      nthreads_(nthreads),
      pixelSizeX_(pixelSizeX),
      pixelSizeY_(pixelSizeY),
      epsilon_(epsilon),
      verbosity_(verbosity) {
  MR_assert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

void WGriddingGridder_Simple::memUsage(size_t &constant,
                                       size_t &per_vis) const {
  // pessimistically assume an oversampling factor of 2
  constant = 4 * width_t_ * height_t_ * sizeof(complex<float>) *
                 2  // sometimes we have two arrays in memory
             + width_t_ * height_t_ * sizeof(float);  // trimmed dirty image
  per_vis =
      sizeof(uint32_t) * 2;  // overestimation, but the best we can do here
}

void WGriddingGridder_Simple::InitializeInversion() {
  img.assign(width_t_ * height_t_, 0);
}

void WGriddingGridder_Simple::AddInversionData(size_t nrows, size_t nchan,
                                               const double *uvw,
                                               const double *freq,
                                               const complex<float> *vis) {
  mav<double, 2> uvw2(uvw, {nrows, 3});
  mav<double, 1> freq2(freq, {nchan});
  mav<complex<float>, 2> ms(vis, {nrows, nchan});
  mav<float, 2> tdirty({width_t_, height_t_});
  mav<float, 2> twgt(nullptr, {0, 0}, false);
  mav<std::uint8_t, 2> tmask(nullptr, {0, 0}, false);
  ms2dirty(uvw2, freq2, ms, twgt, tmask, pixelSizeX_, pixelSizeY_, 0, 0,
           epsilon_, true, nthreads_, tdirty, verbosity_, true);
  for (size_t i = 0; i < width_t_ * height_t_; ++i) img[i] += tdirty[i];
}

void WGriddingGridder_Simple::FinalizeImage(double multiplicationFactor) {
  for (auto &pix : img) pix *= multiplicationFactor;
}

vector<float> WGriddingGridder_Simple::RealImage() {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  vector<float> image(width_ * height_, 0);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      image[(i + dx) + (j + dy) * width_] = img[i * height_t_ + j];
  return image;
}

void WGriddingGridder_Simple::InitializePrediction(vector<float> &&image) {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  img.resize(width_t_ * height_t_);
  MR_assert(image.size() == width_ * height_, "bad image dimensions");
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      img[i * height_t_ + j] = image[(i + dx) + (j + dy) * width_];
}

void WGriddingGridder_Simple::PredictVisibilities(
    size_t nrows, size_t nchan, const double *uvw, const double *freq,
    std::complex<float> *vis) const {
  mav<double, 2> uvw2(uvw, {nrows, 3});
  mav<double, 1> freq2(freq, {nchan});
  mav<complex<float>, 2> ms(vis, {nrows, nchan}, true);
  mav<float, 2> tdirty(img.data(), {width_t_, height_t_});
  mav<float, 2> twgt(nullptr, {0, 0}, false);
  mav<std::uint8_t, 2> tmask(nullptr, {0, 0}, false);
  dirty2ms(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_, pixelSizeY_, 0, 0,
           epsilon_, true, nthreads_, ms, verbosity_, true);
}
