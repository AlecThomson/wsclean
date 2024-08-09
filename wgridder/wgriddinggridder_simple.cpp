#include "wgriddinggridder_simple.h"

#include <complex>
#include <cstddef>
#include <vector>

#include "ducc0/wgridder/wgridder.h"
#include "ducc0/fft/fftnd_impl.h"

#include <LRUCache11.hpp>

#include "../gridding/msgridderbase.h"

using namespace ducc0;

namespace wsclean {

/**
 * VisibilityCallbackBuffer implements a virtual buffer replacement to the
 * `cmav` that would ordinarily be used to pass visibility data into DUCC.
 *
 * Ordinarily the `cmav` that DUCC takes would contain visibilities with facet
 * solutions pre-applied.
 * With VisibilityCallbackBuffer we instead hold in memory a buffer that does
 * not have facet solutions applied.
 * When DUCC requests from the buffer a specific visibility for a specific facet
 * the facet solution is applied "on the fly" and the required value returned.
 * Some internal caching is applied at the row level to help a bit with
 * efficiency.
 */
template <typename TVisibility, GainMode Mode, size_t NPolarizations,
          typename TInfo = ducc0::detail_mav::mav_info<2>>
class VisibilityCallbackBuffer : public TInfo {
 public:
  VisibilityCallbackBuffer(size_t n_rows, size_t n_channels,
                           const aocommon::BandData &selected_band,
                           const std::pair<size_t, size_t> *antenna_buffer,
                           const std::complex<float> *visibility_buffer,
                           const size_t *time_offsets, MSGridderBase *gridder,
                           size_t n_antennas)
      : TInfo({n_rows, n_channels}),
        n_antennas_(n_antennas),
        n_channels_(n_channels),
        n_visibilities_per_row_(n_channels * NPolarizations),
        selected_band_(selected_band),
        antenna_buffer_(antenna_buffer),
        visibility_buffer_(visibility_buffer),
        time_offsets_(time_offsets),
        gridder_(gridder){};

  template <typename Index>
  const TVisibility &raw(Index index) const {
    const size_t row = index / n_channels_;
    const size_t channel = index % n_channels_;
    const std::pair<size_t, size_t> &antenna_pair = antenna_buffer_[row];

    // LRU cache of rows that will grow up to 256 elements in size and then
    // prune itself back to 200 elements
    thread_local static lru11::Cache<size_t,
                                     std::unique_ptr<std::complex<float>[]>>
        visibility_row_cache(200, 56);

    // By caching at a row level we can compute an entire row of modified
    // visibilities at a time instead of individual visibilities as this allows
    // for more efficiency/optimisation in the computation.
    // Subsequent visibility lookups in the same row will hit the cache and not
    // have to recompute; as DUCC always looks up visibilities in a row
    // sequentially this is important for performance and there will be lots of
    // hits.
    // By caching 256 rows we also allow for (less frequent) cache hits when
    // recently requested rows are requested again, this is not as important for
    // performance but still provides some gains.
    if (!visibility_row_cache.contains(row)) {
      std::unique_ptr<std::complex<float>[]> visibility_row(
          new std::complex<float>[n_visibilities_per_row_]);

      std::copy_n(&visibility_buffer_[row * n_visibilities_per_row_],
                  n_visibilities_per_row_, visibility_row.get());

      // We need to pass a cached time_offset into `ApplyCorrections` because
      // it's not capable of calculating  this itself when not called
      // sequentially
      size_t time_offset = time_offsets_[row];
      // We can safely pass null for weight buffer and 0 for time and field_id
      // because these are unused in ModifierBehaviour::kApply mode
      gridder_->ApplyCorrections<Mode, ModifierBehaviour::kApply, false>(
          n_antennas_, visibility_row.get(), selected_band_, nullptr, 0, 0,
          antenna_pair.first, antenna_pair.second, time_offset);

      if constexpr (NPolarizations > 1) {
        internal::CollapseData<NPolarizations>(
            n_channels_, visibility_row.get(), gridder_->Polarization());
      }

      visibility_row_cache.emplace(row, std::move(visibility_row));
    }

    return visibility_row_cache.getRef(row)[channel];
  };
  template <typename... Params>
  const TVisibility operator()(Params... params) const {
    return raw(TInfo::idx(params...));
  }

  // Turn all prefetch operations inside DUCC into null ops
  // As we return by value and are not a persistent buffer prefetching doesn't
  // make sense in this context
  template <typename Index>
  void prefetch_r(Index) const {}
  template <typename Index>
  void prefetch_w(Index) const {}
  template <typename... Params>
  void prefetch_r(Params...) const {}

 private:
  size_t n_antennas_;
  // Number of channels per row of visibilities
  size_t n_channels_;
  size_t n_visibilities_per_row_;
  const aocommon::BandData &selected_band_;
  const std::pair<size_t, size_t> *antenna_buffer_;
  const std::complex<float> *visibility_buffer_;
  /**
   * When applying corrections sequentially a time_offset is calculated by @ref
   * CacheParmResponse() for each row, used for applying the corrections, and
   * then the time_offset for the next row calculated on top of it.
   * As the time_offset is needed when we apply the corrections, and we can't
   * compute it again here without sequentially going through every single row,
   * we have to store all of them in a buffer to be used when we apply the
   * corrections.
   */
  const size_t *time_offsets_;
  MSGridderBase *gridder_;
};

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
template <typename Tms>
inline void WGriddingGridder_Simple<NumT>::AddInversionMs(
    size_t n_rows, const double *uvw, const cmav<double, 1> &freq, Tms &ms) {
  cmav<double, 2> uvw2(uvw, {n_rows, 3});
  vmav<NumT, 2> tdirty({width_t_, height_t_});
  cmav<float, 2> twgt(nullptr, {0, 0});
  cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
  if (!tuning_)
    ms2dirty<NumT, NumT>(uvw2, freq, ms, twgt, tmask, pixelSizeX_, pixelSizeY_,
                         epsilon_, true, nthreads_, tdirty, verbosity_, true,
                         false, sigma_min, sigma_max, -l_shift_, -m_shift_);
  else
    ms2dirty_tuning<NumT, NumT>(uvw2, freq, ms, twgt, tmask, pixelSizeX_,
                                pixelSizeY_, epsilon_, true, nthreads_, tdirty,
                                verbosity_, true, false, sigma_min, sigma_max,
                                -l_shift_, -m_shift_);
  for (size_t i = 0; i < width_t_ * height_t_; ++i) img[i] += tdirty.raw(i);
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
    MSGridderBase *gridder, size_t n_antenna) {
  assert((selected_band.ChannelCount() <= 1) ||
         (frequencies[1] >= frequencies[0]));

  const cmav<double, 1> frequencies2(frequencies,
                                     {selected_band.ChannelCount()});
  // Construct a templated ms:
  //     VisibilityCallbackBuffer<mode, n_polarizations>
  // populated with visibilities and other data and call AddInversionMs(n_rows,
  // uvws, frequencies2, ms) on it.
  AddInversionMs(mode, n_polarizations, n_rows, uvws, frequencies2, n_channels,
                 selected_band, antennas, visibilities, time_offsets, gridder,
                 n_antenna);
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
template <GainMode Mode, typename... Params>
void WGriddingGridder_Simple<NumT>::AddInversionMs(size_t n_polarizations,
                                                   Params... params) {
  switch (n_polarizations) {
    case 0: {
      AddInversionMs<Mode, 0>(params...);
      break;
    }
    case 2: {
      AddInversionMs<Mode, 2>(params...);
      break;
    }
    case 4: {
      AddInversionMs<Mode, 4>(params...);
      break;
    }
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, typename... Params>
void WGriddingGridder_Simple<NumT>::AddInversionMs(size_t n_rows,
                                                   const double *uvw,
                                                   const cmav<double, 1> &freq,
                                                   Params... params) {
  const VisibilityCallbackBuffer<std::complex<float>, Mode, NPolarizations> ms(
      n_rows, params...);
  AddInversionMs(n_rows, uvw, freq, ms);
}

template class WGriddingGridder_Simple<float>;
template class WGriddingGridder_Simple<double>;

}  // namespace wsclean
