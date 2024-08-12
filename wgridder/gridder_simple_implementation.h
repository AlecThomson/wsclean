#ifndef WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_
#define WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_

#include "gridder_simple.h"

#include <complex>
#include <cstddef>
#include <vector>

#include "ducc0/wgridder/wgridder.h"
#include "ducc0/fft/fftnd_impl.h"

#include <LRUCache11.hpp>

#include "../gridding/msgridder.h"

using namespace ducc0;

/*
 * This file contains the implementation of various template methods from @ref
 * WGriddingGridder_Simple They are implemented here instead of directly in the
 * header or a single source file in order to be able to instantiate them in two
 * different source files. This is done because each instantiation must compile
 * the entire class as well as ducc0, which is quite resource intensive even for
 * a single compile. By breaking it into two source files, one for float and one
 * for double we keep time and memory requirements more manageable and allow for
 * better build parrelilisation.
 */

namespace wsclean {

template <typename NumT>
WGriddingGridder_Simple<NumT>::WGriddingGridder_Simple(
    size_t width, size_t height, size_t width_t, size_t height_t,
    double pixelSizeX, double pixelSizeY, double l_shift, double m_shift,
    size_t nthreads, double epsilon, size_t verbosity, bool tuning)
    : width_(width),
      height_(height),
      width_t_(width_t),
      height_t_(height_t),
      nthreads_(nthreads),
      pixelSizeX_(pixelSizeX),
      pixelSizeY_(pixelSizeY),
      l_shift_(l_shift),
      m_shift_(m_shift),
      epsilon_(epsilon),
      verbosity_(verbosity),
      tuning_(tuning) {
  MR_assert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

template <typename NumT>
size_t WGriddingGridder_Simple<NumT>::ConstantMemoryUsage() const {
  // Storage for "grid": pessimistically assume an oversampling factor of 2
  size_t constant = sigma_max * sigma_max * width_t_ * height_t_ *
                    sizeof(std::complex<float>);
  // For prediction, we also need a copy of the dirty image
  constant += width_t_ * height_t_ * sizeof(NumT);  // trimmed dirty image
  return constant;
}

template <typename NumT>
size_t WGriddingGridder_Simple<NumT>::PerVisibilityMemoryUsage() const {
  // Storage for the indexing information is really hard to estimate ...
  // it can go up to 8 bytes per visibility, but this is a really pathological
  // scenario; should typically be below 1 byte/visibility
  return 8;  // Overestimation, but the best we can do here
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::InitializeInversion() {
  img.assign(width_t_ * height_t_, 0);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::AddInversionData(
    size_t n_rows, size_t n_chan, const double *uvw, const double *freq,
    const std::complex<float> *vis) {
  const bool decreasing_freq = (n_chan > 1) && (freq[1] < freq[0]);
  auto freq2(decreasing_freq
                 ? cmav<double, 1>(freq + n_chan - 1, {n_chan}, {-1})
                 : cmav<double, 1>(freq, {n_chan}));
  auto ms(decreasing_freq
              ? cmav<std::complex<float>, 2>(vis + n_chan - 1, {n_rows, n_chan},
                                             {ptrdiff_t(n_chan), -1})
              : cmav<std::complex<float>, 2>(vis, {n_rows, n_chan}));

  AddInversionMs(n_rows, uvw, freq2, ms);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::AddInversionDataWithCorrectionCallback(
    GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
    const double *frequencies, size_t n_channels,
    const aocommon::BandData &selected_band,
    const std::pair<size_t, size_t> *antennas,
    const std::complex<float> *visibilities, const size_t *time_offsets,
    MsGridder *gridder, size_t n_antenna) {
  assert((selected_band.ChannelCount() <= 1) ||
         (frequencies[1] >= frequencies[0]));

  const cmav<double, 1> frequencies2(frequencies,
                                     {selected_band.ChannelCount()});
  // Construct a templated ms:
  //     VisibilityCallbackBuffer<mode, n_polarizations>
  // populated with visibilities and other data and call AddInversionMs(n_rows,
  // uvws, frequencies2, ms) on it.
  AddInversionMs(mode, n_polarizations, n_rows, uvws, std::ref(frequencies2),
                 n_channels, std::ref(selected_band), antennas, visibilities,
                 time_offsets, gridder, n_antenna);
}

template <typename NumT>
template <typename... Params>
void WGriddingGridder_Simple<NumT>::AddInversionMs(GainMode mode,
                                                   Params... params) {
  switch (mode) {
    case GainMode::kXX: {
      AddInversionMs<GainMode::kXX>(params...);
      break;
    }
    case GainMode::kYY: {
      AddInversionMs<GainMode::kYY>(params...);
      break;
    }
    case GainMode::k2VisDiagonal: {
      AddInversionMs<GainMode::k2VisDiagonal>(params...);
      break;
    }
    case GainMode::kTrace: {
      AddInversionMs<GainMode::kTrace>(params...);
      break;
    }
    case GainMode::kFull: {
      AddInversionMs<GainMode::kFull>(params...);
      break;
    }
  }
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::FinalizeImage(double multiplicationFactor) {
  for (auto &pix : img) pix *= multiplicationFactor;
}

template <typename NumT>
std::vector<float> WGriddingGridder_Simple<NumT>::RealImage() {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  std::vector<float> image(width_ * height_,
                           std::numeric_limits<float>::quiet_NaN());
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      image[(i + dx) + (j + dy) * width_] = img[i * height_t_ + j];
  return image;
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::InitializePrediction(
    const float *image_data) {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  img.resize(width_t_ * height_t_);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      img[i * height_t_ + j] = image_data[(i + dx) + (j + dy) * width_];
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::PredictVisibilities(
    size_t n_rows, size_t n_chan, const double *uvw, const double *freq,
    std::complex<float> *vis) const {
  cmav<double, 2> uvw2(uvw, {n_rows, 3});
  bool decreasing_freq = (n_chan > 1) && (freq[1] < freq[0]);
  auto freq2(decreasing_freq
                 ? cmav<double, 1>(freq + n_chan - 1, {n_chan}, {-1})
                 : cmav<double, 1>(freq, {n_chan}));
  auto ms(decreasing_freq
              ? vmav<std::complex<float>, 2>(vis + n_chan - 1, {n_rows, n_chan},
                                             {ptrdiff_t(n_chan), -1})
              : vmav<std::complex<float>, 2>(vis, {n_rows, n_chan}));
  cmav<NumT, 2> tdirty(img.data(), {width_t_, height_t_});
  cmav<float, 2> twgt(nullptr, {0, 0});
  cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
  if (!tuning_)
    dirty2ms<NumT, NumT>(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_,
                         pixelSizeY_, epsilon_, true, nthreads_, ms, verbosity_,
                         true, false, sigma_min, sigma_max, -l_shift_,
                         -m_shift_);
  else
    dirty2ms_tuning<NumT, NumT>(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_,
                                pixelSizeY_, epsilon_, true, nthreads_, ms,
                                verbosity_, true, false, sigma_min, sigma_max,
                                -l_shift_, -m_shift_);
}

}  // namespace wsclean

#endif  //#ifndef WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_
