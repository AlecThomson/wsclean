#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "gridmode.h"

#include <aocommon/banddata.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>
#include <aocommon/imagecoordinates.h>

#include "../structures/observationinfo.h"
#include "../structures/msselection.h"
#include "../structures/weightmode.h"

#include "../msproviders/msreaders/msreader.h"

#include "msprovidercollection.h"
#include "visibilitymodifier.h"
#include "visibilityweightingmode.h"

#include "../main/settings.h"

#include "../scheduling/griddingtaskmanager.h"

#include <aocommon/uvector.h>

#include <mutex>
#include <memory>

namespace schaapcommon {
namespace h5parm {
class H5Parm;
class SolTab;
}  // namespace h5parm
}  // namespace schaapcommon

namespace wsclean {

enum class PsfMode {
  kNone,    // Not a psf, grid the visibilities in the MS
  kSingle,  // Grid generated visibilities for a point source at the centre of
            // the main image
  kDirectionDependent  // Grid generated visibilities for a point source at the
                       // centre of the current facet
};

namespace internal {

template <size_t PolarizationCount>
inline void CollapseData(
    size_t n_channels, std::complex<float>* buffer,
    [[maybe_unused]] aocommon::PolarizationEnum polarization) {
  if constexpr (PolarizationCount == 2) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      buffer[ch] = buffer[ch * PolarizationCount] +
                   buffer[(ch * PolarizationCount + (PolarizationCount - 1))];
    }
  } else if constexpr (PolarizationCount == 4) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      buffer[ch] = aocommon::Polarization::ConvertFromLinear(
          buffer + ch * PolarizationCount, polarization);
    }
  } else
    throw std::runtime_error("Invalid polarization conversion");
}

template <size_t PolarizationCount>
inline void ExpandData(
    size_t n_channels, std::complex<float>* buffer, std::complex<float>* output,
    [[maybe_unused]] aocommon::PolarizationEnum polarization) {
  if constexpr (PolarizationCount == 2) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      output[ch * 2] = buffer[ch];
      output[ch * 2 + 1] = buffer[ch];
    }
  } else if constexpr (PolarizationCount == 4) {
    for (size_t i = 0; i != n_channels; ++i) {
      const size_t ch = n_channels - 1 - i;
      aocommon::Polarization::ConvertToLinear(buffer[ch], polarization,
                                              &output[ch * 4]);
    }
  } else
    throw std::runtime_error("Invalid polarization conversion");
}

}  // namespace internal

class MSGridderBase {
 public:
  MSGridderBase(const Settings& settings,
                MsProviderCollection& ms_provider_collection,
                size_t gridder_index);
  virtual ~MSGridderBase();

  void SetH5Parm(
      const std::vector<schaapcommon::h5parm::H5Parm>& h5parms,
      const std::vector<schaapcommon::h5parm::SolTab*>& first_solutions,
      const std::vector<schaapcommon::h5parm::SolTab*>& second_solutions,
      const std::vector<schaapcommon::h5parm::GainType>& gain_types) {
    visibility_modifier_.SetH5Parm(h5parms, first_solutions, second_solutions,
                                   gain_types);
  }

  /** @return The memory consumption of cached h5 solutions, in bytes. */
  size_t GetCachedH5ParmSize() const {
    return visibility_modifier_.GetCacheParmResponseSize();
  }

  size_t ImageWidth() const { return image_width_; }
  size_t ImageHeight() const { return image_height_; }
  double ImagePadding() const { return image_padding_; }
  double PixelSizeX() const { return settings_.pixelScaleX; }
  double PixelSizeY() const { return settings_.pixelScaleY; }
  size_t ActualWGridSize() const { return actual_w_grid_size_; }

  const Settings& GetSettings() const { return settings_; }

  const std::string& DataColumnName() const { return data_column_name_; }
  bool IsFacet() const { return is_facet_; }
  PsfMode GetPsfMode() const { return psf_mode_; }
  bool DoSubtractModel() const { return do_subtract_model_; }
  bool SmallInversion() const { return small_inversion_; }
  aocommon::PolarizationEnum Polarization() const { return polarization_; }
  WeightMode Weighting() const { return weighting_; }
  const ImageWeights* GetImageWeights() const {
    return precalculated_weight_info_;
  }
  bool IsComplex() const { return is_complex_; }

  VisibilityWeightingMode GetVisibilityWeightingMode() const {
    return visibility_weighting_mode_;
  }
  bool StoreImagingWeights() const { return store_imaging_weights_; }

  void SetFacetIndex(size_t facet_index) { facet_index_ = facet_index; }
  void SetFacetGroupIndex(size_t index) { facet_group_index_ = index; }
  /**
   * @brief In case of facet-based imaging, the model data in the @param
   * MSProvider is reset to zeros in every major cycle, and predicted data
   * should be add-assigned to the model data (_isFacet = true) rather
   * than overwriting it. For standard imaging (_isFacet = false), the model
   * data should be overwritten.
   */
  void SetIsFacet(bool is_facet) { is_facet_ = is_facet; }
  void SetImageWidth(size_t image_width) { image_width_ = image_width; }
  void SetImageHeight(size_t image_height) { image_height_ = image_height; }
  void SetActualWGridSize(size_t actual_w_grid_size) {
    actual_w_grid_size_ = actual_w_grid_size;
  }
  void SetPsfMode(PsfMode psf_mode) { psf_mode_ = psf_mode; }
  void SetImagePadding(double image_padding) { image_padding_ = image_padding; }
  void SetPolarization(aocommon::PolarizationEnum polarization) {
    polarization_ = polarization;
  }
  void SetIsComplex(bool is_complex) { is_complex_ = is_complex; }
  void SetDoSubtractModel(bool do_subtract_model) {
    do_subtract_model_ = do_subtract_model;
  }

  void SetWriterLockManager(GriddingTaskManager* writer_lock_manager) {
    writer_lock_manager_ = writer_lock_manager;
  }

  void SetImageWeights(const ImageWeights* weights) {
    precalculated_weight_info_ = weights;
  }

  /**
   * When processing the first gridder task, the gridder may output more
   * information.
   */
  bool IsFirstTask() const { return is_first_task_; }
  void SetIsFirstTask(bool is_first_task) { is_first_task_ = is_first_task; }

  void SetStoreImagingWeights(bool store_imaging_weights) {
    store_imaging_weights_ = store_imaging_weights;
  }

  virtual void Invert() = 0;

  virtual void Predict(std::vector<aocommon::Image>&& images) = 0;

  virtual std::vector<aocommon::Image> ResultImages() = 0;

  void SetPhaseCentreRA(const double phase_centre_ra) {
    phase_centre_ra_ = phase_centre_ra;
  }
  void SetPhaseCentreDec(const double phase_centre_dec) {
    phase_centre_dec_ = phase_centre_dec;
  }
  double PhaseCentreRA() const { return phase_centre_ra_; }
  double PhaseCentreDec() const { return phase_centre_dec_; }
  void SetLShift(const double l_shift) { l_shift_ = l_shift; }
  void SetMShift(const double m_shift) { m_shift_ = m_shift; }

  void SetMainImageDL(const double main_image_dl) {
    main_image_dl_ = main_image_dl;
  }
  void SetMainImageDM(const double main_image_dm) {
    main_image_dm_ = main_image_dm;
  }

  void SetFacetDirection(double ra, double dec) {
    visibility_modifier_.SetFacetDirection(ra, dec);
  }

  double FacetDirectionRA() const {
    return visibility_modifier_.FacetDirectionRA();
  }
  double FacetDirectionDec() const {
    return visibility_modifier_.FacetDirectionDec();
  }
  double LShift() const { return l_shift_; }
  double MShift() const { return m_shift_; }
  double MainImageDL() const { return main_image_dl_; }
  double MainImageDM() const { return main_image_dm_; }

  /**
   * Deallocate any data that is no longer necessary, but all methods
   * will still return results from the imaging, with the exception of
   * ImageReal/ImageResult().
   */
  virtual void FreeImagingData() {}

  GriddingKernelMode GetGridMode() const { return grid_mode_; }
  void SetGridMode(GriddingKernelMode grid_mode) { grid_mode_ = grid_mode; }

  size_t TrimWidth() const { return trim_width_; }
  size_t TrimHeight() const { return trim_height_; }
  bool HasTrimSize() const { return trim_width_ != 0 || trim_height_ != 0; }
  void SetTrimSize(size_t trim_width, size_t trim_height) {
    trim_width_ = trim_width;
    trim_height_ = trim_height;
  }

  bool HasDenormalPhaseCentre() const {
    return l_shift_ != 0.0 || m_shift_ != 0.0;
  }
  double ImageWeight() const {
    return total_weight_ / GetNVisibilities(gain_mode_);
  }
  /**
   * @return The normalization factor, which is always equal to
   * the image weight in the current implementation.
   *
   * This interface has separate ImageWeight and NormalizationFactor functions
   * since they are conceptually different and the implementation of
   * NormalizationFactor may change in the future.
   */
  double NormalizationFactor() const { return ImageWeight(); }
  double BeamSize() const { return theoretical_beam_size_; }

  /**
   * This is the sum of the weights as given by the measurement set, before the
   * image weighting is applied.
   */
  double VisibilityWeightSum() const { return visibility_weight_sum_; }
  /**
   * The number of visibilities that were gridded.
   */
  size_t GriddedVisibilityCount() const { return gridded_visibility_count_; }
  /**
   * The maximum weight, after having applied the imaging weighting.
   */
  double MaxGriddedWeight() const { return max_gridded_weight_; }
  /**
   * The effective number of visibilities, taking into account imaging weighting
   * and visibility weighting. This number is relative to the "best" visibility:
   * if one visibility with a weight of 10 and 5 visibilities with
   * a weight of 4 were gridded, the effective number of visibilities is
   * (10 + 5 x 4) / 10 = 3
   */
  double EffectiveGriddedVisibilityCount() const {
    return ImageWeight() / MaxGriddedWeight();
  }

  bool HasMetaDataCache() const {
    return (!meta_data_cache_->msDataVector.empty());
  }
  void AllocateMetaDataCache(size_t size) {
    meta_data_cache_->msDataVector.resize(size);
  }
  MetaDataCache::Entry& GetMetaDataCacheItem(size_t index) {
    return meta_data_cache_->msDataVector[index];
  }
  void SetMetaDataCache(std::unique_ptr<MetaDataCache> cache) {
    meta_data_cache_ = std::move(cache);
  }
  std::unique_ptr<MetaDataCache> AcquireMetaDataCache() {
    return std::move(meta_data_cache_);
  }

  void ResetVisibilityModifierCache();

  /**
   * The average squared Mueller correction of all applied corrections.
   * This is the weighted sum of squared Mueller matrices, divided by the sum of
   * weights. It is zero if no corrections are applied.
   * @sa VisibilityModifier::TotalCorrectionSum().
   */
  AverageCorrection GetAverageCorrection() const {
    if (ImageWeight() != 0.0) {
      return visibility_modifier_.TotalCorrectionSum() / ImageWeight();
    } else {
      return AverageCorrection();
    }
  }
  /**
   * The average squared Mueller correction. This is the weighted sum of
   * squared Mueller matrices, divided by the sum of weights.
   * It is zero if not both beam and solution corrections are applied.
   * @sa VisibilityModifier::TotalCorrectionSum().
   */
  AverageCorrection GetAverageBeamCorrection() const {
    if (ImageWeight() != 0.0) {
      return visibility_modifier_.BeamCorrectionSum() / ImageWeight();
    } else {
      return AverageCorrection();
    }
  }

  template <GainMode Mode,
            ModifierBehaviour Behaviour = ModifierBehaviour::kApplyAndSum,
            bool LoadResponse = true>
  void ApplyCorrections(size_t n_antennas, std::complex<float>* visibility_row,
                        const aocommon::BandData& cur_band,
                        float* weight_buffer, double time, size_t field_id,
                        size_t antenna1, size_t antenna2);

  void calculateOverallMetaData(
      const std::vector<MsProviderCollection::FacetData>& ms_facet_data_vector);

  void initializeVisibilityModifierTimes(MsProviderCollection::MsData& msData);

 protected:
  size_t ActualInversionWidth() const { return actual_inversion_width_; }
  size_t ActualInversionHeight() const { return actual_inversion_height_; }
  double ActualPixelSizeX() const { return actual_pixel_size_x_; }
  double ActualPixelSizeY() const { return actual_pixel_size_y_; }

  size_t GetMsCount() const { return ms_count_; }
  MsProviderCollection::MsData& GetMsData(size_t ms_index) {
    return ms_data_vector_[ms_index];
  }
  MsProviderCollection::FacetData& GetMsFacetData(size_t ms_index) {
    return ms_facet_data_vector_[ms_index];
  }

  struct InversionRow {
    double uvw[3];
    std::complex<float>* data;
  };

  /**
   * Initializes MS related data members, i.e. the @c _telescope and the
   * @c _pointResponse data in case a beam is applied on the facets and
   * EveryBeam is available and the @c _predictReader data member in case
   * @c isPredict is true.
   */
  void StartMeasurementSet(const MsProviderCollection::MsData& ms_data,
                           bool is_predict);

  /**
   * Read a row of visibilities from the msprovider, and apply weights, flags
   * and a-terms.
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc).
   *
   * To read the data, this function requires scratch weight and model buffers
   * for storing intermediate values. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   *
   * This function collapses the visibilities in the polarization direction.
   * Gridders that grid a single polarization should use this method instead of
   * @ref GetInstrumentalVisibilities(). The output is stored in the first
   * n_channel elements of the visibility data buffer in @c row_data.
   *
   * @param ms_reader The measurement set provider from which data will be read
   * @param n_antennas The number of antennas
   * @param row_data The resulting weighted data
   * @param cur_band The spectral band currently being imaged
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the full applied weight (i.e. visibility weight * imaging weight).
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   * @param is_selected Per visibility whether that visibility will be gridded
   * in this pass. When the visibility is not gridded, its weight will not be
   * added to the relevant sums (visibility count, weight sum, etc.). This
   * buffer is of size n_chan; i.e. it is not specified per polarization.
   * @param meta_data Metadata that has previously been read from a measurement
   * set provider
   */
  inline void GetCollapsedVisibilities(MSReader& ms_reader, size_t n_antennas,
                                       InversionRow& row_data,
                                       const aocommon::BandData& cur_band,
                                       float* weight_buffer,
                                       std::complex<float>* model_buffer,
                                       const bool* is_selected,
                                       const MSProvider::MetaData& meta_data) {
    ReadVisibilities(ms_reader, row_data, weight_buffer, model_buffer);

    CollapseVisibilities(n_antennas, row_data, cur_band, weight_buffer,
                         model_buffer, is_selected, meta_data);

    if (StoreImagingWeights())
      ms_reader.WriteImagingWeights(scratch_image_weights_.data());
  }

  /**
   * Same as @ref GetCollapsedVisibilities(), but without collapsing the
   * polarization direction. This implies that the output visibility buffer in
   * the row_data structure will contain n_channel x n_polarization elements.
   */
  template <size_t PolarizationCount>
  inline void GetInstrumentalVisibilities(
      MSReader& ms_reader, size_t n_antennas, InversionRow& row_data,
      const aocommon::BandData& cur_band, float* weight_buffer,
      std::complex<float>* model_buffer, const bool* is_selected,
      const MSProvider::MetaData& meta_data) {
    ReadVisibilities(ms_reader, row_data, weight_buffer, model_buffer);

    CalculateWeights<PolarizationCount>(row_data, cur_band, weight_buffer,
                                        model_buffer, is_selected);

    ApplyWeightsAndCorrections(n_antennas, row_data, cur_band, weight_buffer,
                               meta_data);

    if (StoreImagingWeights())
      ms_reader.WriteImagingWeights(scratch_image_weights_.data());
  }

  /**
   * @brief Apply corrections as well as visibility and imaging weights
   * Also computes the weight corresponding to the
   * combined effect of the corrections.
   *
   * Requires `scratch_image_weights_` to be populated which is usually done by
   * calling @ref CalculateWeights()
   */
  void ApplyWeightsAndCorrections(size_t n_antennas, InversionRow& row_data,
                                  const aocommon::BandData& cur_band,
                                  float* weight_buffer,
                                  const MSProvider::MetaData& meta_data);

  /**
   * Read a row of visibility and weights from the msprovider
   *
   * Use this function to correctly populate an InversionRow structure and an
   * accompanying weight_buffer and model_buffer before calling @ref
   * CollapseVisibilities() or @ref ApplyWeightsAndCorrections()
   *
   * @param ms_reader The measurement set provider from which data will be read
   * @param row_data The caller must set this object up to point at the desired
   * portion of an allocated buffer into which the visibilities will be read.
   * After returning from this call the uvw paramater of this object will be
   * populated with the (u/v/w)InM values of `meta_data`
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the weights from `ms_reader`
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   */
  inline void ReadVisibilities(MSReader& ms_reader, InversionRow& row_data,
                               float* weight_buffer,
                               std::complex<float>* model_buffer) {
    if (GetPsfMode() == PsfMode::kNone) {
      ms_reader.ReadData(row_data.data);
    }
    if (DoSubtractModel()) {
      ms_reader.ReadModel(model_buffer);
    }
    ms_reader.ReadWeights(weight_buffer);
  }

  /**
   * Apply weights, flags and a-terms to a row of visibility data that has been
   * read by @ref ReadVisibilities() and collapse in the polarization direction
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc).
   *
   * To read the data, this function requires scratch weight and model buffers
   * for storing intermediate values. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   *
   * This function collapses the visibilities in the polarization direction.
   * Gridders that grid a single polarization should use this method instead of
   * @ref ApplyWeightsAndCorrections(). The output is stored in the first
   * n_channel elements of the visibility data buffer in @c row_data.
   *
   * Normally set to one when imaging a single
   * polarization. It may be set to 2 or 4 for IDG as it images multiple
   * polarizations at once, and it may be set to 2 or 4 when applying
   * solutions.
   * @param row_data The resulting weighted data
   * @param cur_band The spectral band currently being imaged
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the full applied weight (i.e. visibility weight * imaging weight).
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   * @param is_selected Per visibility whether that visibility will be gridded
   * in this pass. When the visibility is not gridded, its weight will not be
   * added to the relevant sums (visibility count, weight sum, etc.). This
   * buffer is of size n_chan; i.e. it is not specified per polarization.
   * @param meta_data Metadata that has previously been read from a measurement
   * set provider
   */
  inline void CollapseVisibilities(size_t n_antennas, InversionRow& row_data,
                                   const aocommon::BandData& cur_band,
                                   float* weight_buffer,
                                   std::complex<float>* model_buffer,
                                   const bool* is_selected,
                                   const MSProvider::MetaData& meta_data) {
    switch (n_vis_polarizations_) {
      case 1:
        CalculateWeights<1>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, cur_band,
                                   weight_buffer, meta_data);
        break;
      case 2:
        CalculateWeights<2>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, cur_band,
                                   weight_buffer, meta_data);
        internal::CollapseData<2>(cur_band.ChannelCount(), row_data.data,
                                  Polarization());
        break;
      case 4:
        CalculateWeights<4>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, cur_band,
                                   weight_buffer, meta_data);
        internal::CollapseData<4>(cur_band.ChannelCount(), row_data.data,
                                  Polarization());
        break;
    }
  }

  /**
   * Write (modelled) visibilities to MS, provides an interface to
   * MSProvider::WriteModel(). Any active facet beam or solution corrections
   * are applied. Method is templated on the number of
   * polarizations (1, 2 or 4). The gain_mode can be used to
   * select an entry or entries from the gain matrix that should be used for the
   * correction.
   * @param buffer should on entry contain n_channels visibilities to be
   * written.
   */
  void WriteCollapsedVisibilities(MSProvider& ms_provider, size_t n_antennas,
                                  const aocommon::BandData& cur_band,
                                  std::complex<float>* buffer,
                                  MSProvider::MetaData& meta_data) {
    switch (n_vis_polarizations_) {
      case 1:
        WriteInstrumentalVisibilities(ms_provider, n_antennas, cur_band, buffer,
                                      meta_data);
        break;
      case 2:
        internal::ExpandData<2>(cur_band.ChannelCount(), buffer,
                                scratch_model_data_.data(), Polarization());
        WriteInstrumentalVisibilities(ms_provider, n_antennas, cur_band,
                                      scratch_model_data_.data(), meta_data);
        break;
      case 4:
        internal::ExpandData<4>(cur_band.ChannelCount(), buffer,
                                scratch_model_data_.data(), Polarization());
        WriteInstrumentalVisibilities(ms_provider, n_antennas, cur_band,
                                      scratch_model_data_.data(), meta_data);
        break;
    }
  }

  /**
   * Similar to @ref WriteCollapsedVisibilities(), but assumes the input are
   * instrumental visibilities.
   * @param buffer n_polarizations x n_channels entries, which are the
   * instrumental visibilities.
   */
  void WriteInstrumentalVisibilities(MSProvider& ms_rovider, size_t n_antennas,
                                     const aocommon::BandData& cur_band,
                                     std::complex<float>* buffer,
                                     MSProvider::MetaData& meta_data);

  virtual size_t getSuggestedWGridSize() const = 0;

  void resetVisibilityCounters() {
    gridded_visibility_count_ = 0;
    total_weight_ = 0.0;
    max_gridded_weight_ = 0.0;
    visibility_weight_sum_ = 0.0;
  }

  template <size_t PolarizationCount>
  static void rotateVisibilities(const aocommon::BandData& bandData,
                                 double shiftFactor,
                                 std::complex<float>* dataIter);

  void ReadPredictMetaData(MSProvider::MetaData& meta_data);

  /**
   * The largest w value present in the data (after applying any selections),
   * in units of number of wavelengths. It is initialized by
   * initializeMSDataVector() and is undefined beforehand.
   */
  double MaximumW() const { return max_w_; }
  double MinimumW() const { return min_w_; }

 private:
  bool hasWGridSize() const { return w_grid_size_ != 0; }
  void initializePointResponse(const MsProviderCollection::MsData& msData);

  template <size_t PolarizationCount>
  void CalculateWeights(InversionRow& row_data,
                        const aocommon::BandData& cur_band,
                        float* weight_buffer, std::complex<float>* model_buffer,
                        const bool* is_selected);

  /**
   * @brief Applies the selected visibility modifier (selected by Mode)
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @tparam Behaviour See @ref ModifierBehaviour and @ref @ref
   * ApplyConjugatedParmResponse for more information
   * @tparam LoadResponse This should always be true unless the calling code
   * knows the response has already been loaded previously, e.g. if we first
   * call `ApplyCorrections<Mode, kSum, true>(...)` we can then afterwards call
   * `ApplyCorrection<Mode, kApply, false>(...)` for the same values
   */
  template <GainMode Mode,
            ModifierBehaviour Behaviour = ModifierBehaviour::kApplyAndSum,
            bool LoadResponse = true>
  void ApplyCorrections(size_t n_antennas, std::complex<float>* visibility_row,
                        const aocommon::BandData& cur_band,
                        float* weight_buffer,
                        const MSProvider::MetaData& meta_data) {
    ApplyCorrections<Mode, Behaviour, LoadResponse>(
        n_antennas, visibility_row, cur_band, weight_buffer, meta_data.time,
        meta_data.fieldId, meta_data.antenna1, meta_data.antenna2);
  }

  /**
   * @brief Apply visibility and imaging weights
   * Requires `scratch_image_weights_` to be populated which is usually done by
   * calling @ref CalculateWeights()
   */
  template <GainMode Mode>
  void ApplyWeights(std::complex<float>* visibility_row,
                    const size_t channel_count, float* weight_buffer);

  template <GainMode Mode>
  void WriteInstrumentalVisibilities(MSProvider& ms_provider, size_t n_antennas,
                                     const aocommon::BandData& cur_band,
                                     std::complex<float>* buffer,
                                     MSProvider::MetaData& meta_data);

  /**
   * @brief Applies both the conjugated h5 parm
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedH5Parm(MSReader& ms_reader,
                             const std::vector<std::string>& antenna_names,
                             InversionRow& row_data,
                             const aocommon::BandData& cur_band,
                             const float* weight_buffer,
                             bool apply_forward = false);

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Applies the conjugated facet beam to the visibilities and computes
   * the weight corresponding to the combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */

  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedFacetBeam(MSReader& ms_reader, InversionRow& row_data,
                                const aocommon::BandData& cur_band,
                                const float* weight_buffer,
                                bool apply_forward = false);

  /**
   * @brief Applies both the conjugated facet beam and the conjugated h5 parm
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedFacetDdEffects(
      MSReader& ms_reader, const std::vector<std::string>& antenna_names,
      InversionRow& row_data, const aocommon::BandData& cur_band,
      const float* weight_buffer, bool apply_forward = false);
#endif  // HAVE_EVERYBEAM

  // Contains all MS metadata as well as MS specific gridding data
  size_t ms_count_;
  std::vector<MsProviderCollection::MsData>& ms_data_vector_;
  std::vector<MsProviderCollection::FacetData>& ms_facet_data_vector_;

  const Settings& settings_;
  std::unique_ptr<MetaDataCache> meta_data_cache_;
  size_t actual_inversion_width_ = 0;
  size_t actual_inversion_height_ = 0;
  double actual_pixel_size_x_ = 0.0;
  double actual_pixel_size_y_ = 0.0;
  double phase_centre_ra_ = 0.0;
  double phase_centre_dec_ = 0.0;
  double l_shift_ = 0.0;
  double m_shift_ = 0.0;
  double main_image_dl_ = 0.0;
  double main_image_dm_ = 0.0;
  size_t facet_index_ = 0;
  /// @p _facetGroupIndex and @p _msIndex in conjunction with the @p
  /// GetMsCount() determine the index in the _writerGroupLocks
  /// vector, having size FacetGroupCount() * GetMsCount(). These
  /// variable are only relevant for prediction.
  size_t facet_group_index_ = 0;
  size_t original_ms_index_ = 0;
  /// @see SetIsFacet()
  bool is_facet_ = false;
  double image_padding_ = 1.0;
  size_t image_width_ = 0.0;
  size_t image_height_ = 0.0;
  size_t trim_width_ = 0;
  size_t trim_height_ = 0;
  size_t w_grid_size_ = 0;
  size_t actual_w_grid_size_ = 0;
  std::string data_column_name_;
  PsfMode psf_mode_ = PsfMode::kNone;
  bool do_subtract_model_ = false;
  bool small_inversion_ = true;
  double max_w_ = 0.0;
  double min_w_ = 0.0;
  /// A fractional value that, when non-zero, places a limit on the w-value of
  /// gridded visibilities. Visibilities outside the limit are skipped.
  double w_limit_ = 0.0;
  const ImageWeights* precalculated_weight_info_ = nullptr;
  aocommon::PolarizationEnum polarization_ = aocommon::Polarization::StokesI;
  size_t n_vis_polarizations_ = 1;
  GainMode gain_mode_ = GainMode::kTrace;
  bool is_complex_ = false;
  WeightMode weighting_ = WeightMode(WeightMode::UniformWeighted);
  bool is_first_task_ = false;
  VisibilityWeightingMode visibility_weighting_mode_ =
      VisibilityWeightingMode::NormalVisibilityWeighting;
  GriddingKernelMode grid_mode_ = GriddingKernelMode::KaiserBessel;
  bool store_imaging_weights_ = false;
  double theoretical_beam_size_ = 0.0;

  size_t gridded_visibility_count_ = 0;
  double total_weight_ = 0.0;
  double max_gridded_weight_ = 0.0;
  double visibility_weight_sum_ = 0.0;

  aocommon::UVector<float> scratch_image_weights_;
  /// Initialized in StartMeasurementSet(), used in WriteCollapsedVisibilities()
  /// to expand visibilities into.
  aocommon::UVector<std::complex<float>> scratch_model_data_;

  std::unique_ptr<MSReader> predict_reader_;
  GriddingTaskManager* writer_lock_manager_ = nullptr;

  VisibilityModifier visibility_modifier_;
};

template <GainMode Mode, ModifierBehaviour Behaviour, bool LoadResponse>
inline void MSGridderBase::ApplyCorrections(size_t n_antennas,
                                            std::complex<float>* visibility_row,
                                            const aocommon::BandData& cur_band,
                                            float* weight_buffer, double time,
                                            size_t field_id, size_t antenna1,
                                            size_t antenna2) {
  if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    const bool apply_forward = GetPsfMode() == PsfMode::kDirectionDependent;
    if (apply_beam && visibility_modifier_.HasH5Parm()) {
#ifdef HAVE_EVERYBEAM
      // Load and apply (in conjugate) both the beam and the h5parm solutions
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheBeamResponse(time, field_id, cur_band);
        visibility_modifier_.CacheParmResponse(time, cur_band,
                                               original_ms_index_);
      }
      visibility_modifier_.ApplyConjugatedDual<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights_.data(),
          cur_band.ChannelCount(), n_antennas, antenna1, antenna2,
          original_ms_index_, apply_forward);
    } else if (apply_beam) {
      // Load and apply only the conjugate beam
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheBeamResponse(time, field_id, cur_band);
      }
      visibility_modifier_.ApplyConjugatedBeamResponse<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights_.data(),
          cur_band.ChannelCount(), antenna1, antenna2, apply_forward);

#endif  // HAVE_EVERYBEAM
    } else if (visibility_modifier_.HasH5Parm()) {
      // Load and apply the h5parm solutions
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheParmResponse(time, cur_band,
                                               original_ms_index_);
      }
      visibility_modifier_.ApplyConjugatedParmResponse<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights_.data(),
          original_ms_index_, cur_band.ChannelCount(), n_antennas, antenna1,
          antenna2, apply_forward);
    }
  }
}

template <GainMode Mode>
inline void MSGridderBase::ApplyWeights(std::complex<float>* visibility_row,
                                        const size_t channel_count,
                                        float* weight_buffer) {
  const size_t n_pols = GetNVisibilities(Mode);

  for (size_t channel = 0; channel < channel_count; channel++) {
    for (size_t pol = 0; pol < n_pols; pol++) {
      size_t i = channel * n_pols + pol;

      const float cumWeight =
          weight_buffer[i] * scratch_image_weights_[channel];
      // We can use the boolean for computation instead of an if-condition
      // within the loop. This allows the inner part of the loop to be
      // autovectorized more easily.
      const bool noWeight = cumWeight == 0.0;
      if (pol == 0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        visibility_weight_sum_ += weight_buffer[i] * noWeight;
        max_gridded_weight_ =
            std::max(static_cast<double>(cumWeight), max_gridded_weight_);
        gridded_visibility_count_ += !noWeight;
      }
      // Total weight includes imaging weights
      total_weight_ += cumWeight;
      weight_buffer[i] = cumWeight;
      visibility_row[i] *= cumWeight;
    }
  }
}

// Apply corrections as well as visibility and imaging weights
inline void MSGridderBase::ApplyWeightsAndCorrections(
    size_t n_antennas, InversionRow& row_data,
    const aocommon::BandData& cur_band, float* weight_buffer,
    const MSProvider::MetaData& meta_data) {
  switch (gain_mode_) {
    case GainMode::kXX:
      ApplyCorrections<GainMode::kXX>(n_antennas, row_data.data, cur_band,
                                      weight_buffer, meta_data);
      ApplyWeights<GainMode::kXX>(row_data.data, cur_band.ChannelCount(),
                                  weight_buffer);
      break;
    case GainMode::kYY:
      ApplyCorrections<GainMode::kYY>(n_antennas, row_data.data, cur_band,
                                      weight_buffer, meta_data);
      ApplyWeights<GainMode::kYY>(row_data.data, cur_band.ChannelCount(),
                                  weight_buffer);
      break;
    case GainMode::kTrace:
      ApplyCorrections<GainMode::kTrace>(n_antennas, row_data.data, cur_band,
                                         weight_buffer, meta_data);
      ApplyWeights<GainMode::kTrace>(row_data.data, cur_band.ChannelCount(),
                                     weight_buffer);
      break;
    case GainMode::k2VisDiagonal:
      ApplyCorrections<GainMode::k2VisDiagonal>(
          n_antennas, row_data.data, cur_band, weight_buffer, meta_data);
      ApplyWeights<GainMode::k2VisDiagonal>(
          row_data.data, cur_band.ChannelCount(), weight_buffer);
      break;
    case GainMode::kFull:
      ApplyCorrections<GainMode::kFull>(n_antennas, row_data.data, cur_band,
                                        weight_buffer, meta_data);
      ApplyWeights<GainMode::kFull>(row_data.data, cur_band.ChannelCount(),
                                    weight_buffer);
      break;
    default:
      throw std::runtime_error(
          "Invalid combination of visibility polarizations and gain mode");
  }
}

inline void MSGridderBase::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& cur_band, std::complex<float>* buffer,
    MSProvider::MetaData& meta_data) {
  switch (gain_mode_) {
    case GainMode::kXX:
      WriteInstrumentalVisibilities<GainMode::kXX>(ms_provider, n_antennas,
                                                   cur_band, buffer, meta_data);
      break;
    case GainMode::kYY:
      WriteInstrumentalVisibilities<GainMode::kYY>(ms_provider, n_antennas,
                                                   cur_band, buffer, meta_data);
      break;
    case GainMode::kTrace:
      WriteInstrumentalVisibilities<GainMode::kTrace>(
          ms_provider, n_antennas, cur_band, buffer, meta_data);
      break;
    case GainMode::k2VisDiagonal:
      WriteInstrumentalVisibilities<GainMode::k2VisDiagonal>(
          ms_provider, n_antennas, cur_band, buffer, meta_data);
      break;
    case GainMode::kFull:
      WriteInstrumentalVisibilities<GainMode::kFull>(
          ms_provider, n_antennas, cur_band, buffer, meta_data);
      break;
  }
}

}  // namespace wsclean

#endif
