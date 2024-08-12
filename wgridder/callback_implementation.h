#ifndef WSCLEAN_WGRIDDER_CALLBACK_IMPL_H_
#define WSCLEAN_WGRIDDER_CALLBACK_IMPL_H_

/*
 * This file contains the implementation of @ref VisibilityCallbackBuffer and
 * some of the methods from @ref WGriddingGridder_Simple used to
 * instantiate/call @ref VisibilityCallbackBuffer. They are implemented here
 * instead of directly in wgriddingsimple.h or in a source file so that they can
 * be instantiated in multiple different source files each with a different
 * explicit instantiation. This is because each instantiation must compile
 * ducc0. Which is quite resource intensive even for a single compile, and as we
 * have `NumT * GainMode * Polarization` instantiations or `2 * 5 * 3 = 30`
 * instantiations. We instantiate in 10 different files with 3 instantiations
 * per file which keeps time and memory requirements are kept more manageable
 */

class MsGridder;

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
                           const size_t *time_offsets, MsGridder *gridder,
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
  MsGridder *gridder_;
};

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, typename... Params>
void WGriddingGridder_Simple<NumT>::AddInversionMs(
    size_t n_rows, const double *uvw, const ducc0::cmav<double, 1> &freq,
    Params... params) {
  const VisibilityCallbackBuffer<std::complex<float>, Mode, NPolarizations> ms(
      n_rows, params...);
  AddInversionMs(n_rows, uvw, freq, ms);
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

#endif  // WSCLEAN_WGRIDDER_CALLBACK_IMPL_H_
