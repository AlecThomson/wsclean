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

#include "visibilitymodifier.h"
#include "visibilityweightingmode.h"

#include "../main/settings.h"

#include "../scheduling/griddingtaskmanager.h"

#include <aocommon/uvector.h>

#include <map>
#include <mutex>
#include <memory>

namespace schaapcommon {
namespace h5parm {
class H5Parm;
class SolTab;
}  // namespace h5parm
}  // namespace schaapcommon

enum class PsfMode {
  kNone,    // Not a psf, grid the visibilities in the MS
  kSingle,  // Grid generated visibilities for a point source at the centre of
            // the main image
  kDirectionDependent  // Grid generated visibilities for a point source at the
                       // centre of the current facet
};

namespace internal {

template <size_t PolarizationCount>
inline void CollapseData(size_t n_channels, std::complex<float>* buffer) {
  // PolarizationCount is either 2 or 4. In the case it is two, we need to add
  // element 0 and 1. In the case it is four, we need to add element 0 and 3.
  for (size_t ch = 0; ch != n_channels; ++ch) {
    buffer[ch] = buffer[ch * PolarizationCount] +
                 buffer[(ch * PolarizationCount + (PolarizationCount - 1))];
  }
}

template <size_t PolarizationCount>
inline void ExpandData(size_t n_channels, std::complex<float>* buffer,
                       std::complex<float>* output);

template <>
inline void ExpandData<2>(size_t n_channels, std::complex<float>* buffer,
                          std::complex<float>* output) {
  for (size_t ch = 0; ch != n_channels; ++ch) {
    output[ch * 2] = buffer[ch];
    output[ch * 2 + 1] = buffer[ch];
  }
}

template <>
inline void ExpandData<4>(size_t n_channels, std::complex<float>* buffer,
                          std::complex<float>* output) {
  for (size_t i = 0; i != n_channels; ++i) {
    const size_t ch = n_channels - 1 - i;
    output[ch * 4] = buffer[ch];
    output[ch * 4 + 1] = 0.0;
    output[ch * 4 + 2] = 0.0;
    output[ch * 4 + 3] = buffer[ch];
  }
}

}  // namespace internal

class MSGridderBase {
 public:
  MSGridderBase(const Settings& settings);
  virtual ~MSGridderBase();

  size_t ImageWidth() const { return image_width_; }
  size_t ImageHeight() const { return image_height_; }
  double ImagePadding() const { return image_padding_; }
  double PixelSizeX() const { return settings_.pixelScaleX; }
  double PixelSizeY() const { return settings_.pixelScaleY; }
  size_t ActualWGridSize() const { return actual_w_grid_size_; }

  void ClearMeasurementSetList() {
    measurement_sets_.clear();
    selections_.clear();
  }
  const Settings& GetSettings() const { return settings_; }
  MSProvider& MeasurementSet(size_t index) const {
    return *measurement_sets_[index];
  }
  const MSSelection& Selection(size_t index) const {
    return selections_[index];
  }
  size_t MeasurementSetCount() const { return measurement_sets_.size(); }
  void AddMeasurementSet(std::unique_ptr<MSProvider> ms_provider,
                         const MSSelection& selection) {
    measurement_sets_.push_back(std::move(ms_provider));
    selections_.push_back(selection);
  }

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

  double StartTime() const { return start_time_; }
  bool HasDenormalPhaseCentre() const {
    return l_shift_ != 0.0 || m_shift_ != 0.0;
  }
  double ImageWeight() const { return total_weight_; }
  /**
   * @return The normalization factor, which is always equal to
   * the image weight in the current implementation.
   *
   * This interface has separate ImageWeight and NormalizationFactor functions
   * since they are conceptually different and the implementation of
   * NormalizationFactor may change in the future.
   */
  double NormalizationFactor() const { return total_weight_; }
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
    return totalWeight() / MaxGriddedWeight();
  }

  void SetMetaDataCache(std::unique_ptr<MetaDataCache> cache) {
    meta_data_cache_ = std::move(cache);
  }
  std::unique_ptr<MetaDataCache> AcquireMetaDataCache() {
    return std::move(meta_data_cache_);
  }

  double AverageCorrection() const {
    return visibility_modifier_.CorrectionSum() / total_weight_;
  }
  double AverageH5Correction() const {
    return visibility_modifier_.H5CorrectionSum() / total_weight_;
  }

 protected:
  size_t ActualInversionWidth() const { return actual_inversion_width_; }
  size_t ActualInversionHeight() const { return actual_inversion_height_; }
  double ActualPixelSizeX() const { return actual_pixel_size_x_; }
  double ActualPixelSizeY() const { return actual_pixel_size_y_; }

  struct MSData {
   public:
    MSProvider* ms_provider = nullptr;
    size_t msIndex = 0;
    size_t dataDescId = 0;
    aocommon::BandData bandData;
    size_t startChannel = 0;
    size_t endChannel = 0;
    size_t matchingRows = 0;
    size_t totalRowsProcessed = 0;
    double minW = 0.0;
    double maxW = 0.0;
    double maxWWithFlags = 0.0;
    double maxBaselineUVW = 0.0;
    double maxBaselineInM = 0.0;
    size_t rowStart = 0;
    size_t rowEnd = 0;
    double integrationTime = 0.0;
    std::vector<std::string> antenna_names;

    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(bandData, startChannel, endChannel);
    }
  };

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
  void StartMeasurementSet(const MSGridderBase::MSData& ms_data,
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
   * @param antenna_names The antenna names
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
  inline void GetCollapsedVisibilities(
      MSReader& ms_reader, const std::vector<std::string>& antenna_names,
      InversionRow& row_data, const aocommon::BandData& cur_band,
      float* weight_buffer, std::complex<float>* model_buffer,
      const bool* is_selected, const MSProvider::MetaData& meta_data) {
    ReadVisibilities(ms_reader, row_data, weight_buffer, model_buffer);

    CollapseVisibilities(antenna_names, row_data, cur_band, weight_buffer,
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
      MSReader& ms_reader, const std::vector<std::string>& antenna_names,
      InversionRow& row_data, const aocommon::BandData& cur_band,
      float* weight_buffer, std::complex<float>* model_buffer,
      const bool* is_selected, const MSProvider::MetaData& meta_data) {
    ReadVisibilities(ms_reader, row_data, weight_buffer, model_buffer);

    CalculateWeights<PolarizationCount>(row_data, cur_band, weight_buffer,
                                        model_buffer, is_selected);

    ApplyWeightsAndCorrections<PolarizationCount>(
        antenna_names, row_data, cur_band, weight_buffer, meta_data);

    if (StoreImagingWeights())
      ms_reader.WriteImagingWeights(scratch_image_weights_.data());
  }

  inline void CalculateWeights(double* uvw_buffer,
                               std::complex<float>* visibility_buffer,
                               const aocommon::BandData& cur_band,
                               float* weight_buffer,
                               std::complex<float>* model_buffer,
                               const bool* is_selected) {
    switch (n_vis_polarizations_) {
      case 1:
        CalculateWeights<1>(uvw_buffer, visibility_buffer, cur_band,
                            weight_buffer, model_buffer, is_selected);
        break;
      case 2:
        CalculateWeights<2>(uvw_buffer, visibility_buffer, cur_band,
                            weight_buffer, model_buffer, is_selected);
        break;
      case 4:
        CalculateWeights<4>(uvw_buffer, visibility_buffer, cur_band,
                            weight_buffer, model_buffer, is_selected);
        break;
    }
  }

  template <size_t PolarizationCount>
  std::complex<float>& InlineApplyWeightsAndCorrections(
      const std::vector<std::string>& antenna_names,
      const aocommon::BandData& cur_band, float* scratch_weight_buffer,
      float* weight_buffer, std::complex<float>* visibility_buffer,
      size_t visibility_index, size_t row, size_t channel,
      const MSProvider::MetaData& meta_data);

  template <size_t PolarizationCount>
  void ApplyWeightsAndCorrections(const std::vector<std::string>& antenna_names,
                                  InversionRow& row_data,
                                  const aocommon::BandData& cur_band,
                                  float* weight_buffer,
                                  const MSProvider::MetaData& meta_data);

  void PreCacheCorrections(const MSProvider::MetaData& metaData,
                           const std::vector<std::string>& antenna_names,
                           const aocommon::BandData& curBand);

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
    ReadVisibilities(ms_reader, row_data.data, weight_buffer, model_buffer);
  }
  inline void ReadVisibilities(MSReader& ms_reader,
                               std::complex<float>* visibility_buffer,
                               float* weight_buffer,
                               std::complex<float>* model_buffer) {
    if (GetPsfMode() == PsfMode::kNone) {
      ms_reader.ReadData(visibility_buffer);
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
   * solulutions.
   * @param antenna_names The antenna names
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
  inline void CollapseVisibilities(
      const std::vector<std::string>& antenna_names, InversionRow& row_data,
      const aocommon::BandData& cur_band, float* weight_buffer,
      std::complex<float>* model_buffer, const bool* is_selected,
      const MSProvider::MetaData& meta_data) {
    switch (n_vis_polarizations_) {
      case 1:
        CalculateWeights<1>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections<1>(antenna_names, row_data, cur_band,
                                      weight_buffer, meta_data);
        break;
      case 2:
        CalculateWeights<2>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections<2>(antenna_names, row_data, cur_band,
                                      weight_buffer, meta_data);
        internal::CollapseData<2>(cur_band.ChannelCount(), row_data.data);
        break;
      case 4:
        CalculateWeights<4>(row_data, cur_band, weight_buffer, model_buffer,
                            is_selected);
        ApplyWeightsAndCorrections<4>(antenna_names, row_data, cur_band,
                                      weight_buffer, meta_data);
        internal::CollapseData<4>(cur_band.ChannelCount(), row_data.data);
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
  void WriteCollapsedVisibilities(MSProvider& ms_provider,
                                  const std::vector<std::string>& antenna_names,
                                  const aocommon::BandData& cur_band,
                                  std::complex<float>* buffer,
                                  MSProvider::MetaData& meta_data) {
    switch (n_vis_polarizations_) {
      case 1:
        WriteInstrumentalVisibilities<1>(ms_provider, antenna_names, cur_band,
                                         buffer, meta_data);
        break;
      case 2:
        internal::ExpandData<2>(cur_band.ChannelCount(), buffer,
                                scratch_model_data_.data());
        WriteInstrumentalVisibilities<2>(ms_provider, antenna_names, cur_band,
                                         scratch_model_data_.data(), meta_data);
        break;
      case 4:
        internal::ExpandData<4>(cur_band.ChannelCount(), buffer,
                                scratch_model_data_.data());
        WriteInstrumentalVisibilities<4>(ms_provider, antenna_names, cur_band,
                                         scratch_model_data_.data(), meta_data);
        break;
    }
  }

  /**
   * Similar to @ref WriteCollapsedVisibilities(), but assumes the input are
   * instrumental visibilities.
   * @param Buffer with PolarizationCount x n_channels entries, which are the
   * instrumental visibilities.
   */
  template <size_t PolarizationCount>
  void WriteInstrumentalVisibilities(
      MSProvider& ms_rovider, const std::vector<std::string>& antenna_names,
      const aocommon::BandData& cur_band, std::complex<float>* buffer,
      MSProvider::MetaData& meta_data);

  virtual size_t getSuggestedWGridSize() const = 0;

  void resetVisibilityCounters() {
    gridded_visibility_count_ = 0;
    total_weight_ = 0.0;
    max_gridded_weight_ = 0.0;
    visibility_weight_sum_ = 0.0;
  }

  double totalWeight() const { return total_weight_; }

  void initializeMSDataVector(std::vector<MSData>& msDataVector);

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
  static std::vector<std::string> getAntennaNames(
      const casacore::MSAntenna& msAntenna);

  void resetMetaData() { has_frequencies_ = false; }

  void calculateMSLimits(const aocommon::BandData& selectedBand,
                         double startTime) {
    if (has_frequencies_) {
      freq_low_ = std::min(freq_low_, selectedBand.LowestFrequency());
      freq_high_ = std::max(freq_high_, selectedBand.HighestFrequency());
      band_start_ = std::min(band_start_, selectedBand.BandStart());
      band_end_ = std::max(band_end_, selectedBand.BandEnd());
      start_time_ = std::min(start_time_, startTime);
    } else {
      freq_low_ = selectedBand.LowestFrequency();
      freq_high_ = selectedBand.HighestFrequency();
      band_start_ = selectedBand.BandStart();
      band_end_ = selectedBand.BandEnd();
      start_time_ = startTime;
      has_frequencies_ = true;
    }
  }

  template <size_t NPolInMSProvider>
  void calculateWLimits(MSGridderBase::MSData& msData);

  void initializeMeasurementSet(MSGridderBase::MSData& msData,
                                MetaDataCache::Entry& cacheEntry,
                                bool isCacheInitialized);

  void calculateOverallMetaData(const std::vector<MSData>& msDataVector);
  bool hasWGridSize() const { return w_grid_size_ != 0; }
  void initializeBandData(const casacore::MeasurementSet& ms, MSData& msData);
  void initializePointResponse(const MSData& msData);

  template <size_t PolarizationCount>
  inline void CalculateWeights(InversionRow& rowData,
                               const aocommon::BandData& curBand,
                               float* weightBuffer,
                               std::complex<float>* modelBuffer,
                               const bool* isSelected) {
    CalculateWeights<PolarizationCount>(rowData.uvw, rowData.data, curBand,
                                        weightBuffer, modelBuffer, isSelected);
  }

  template <size_t PolarizationCount>
  void CalculateWeights(double* uvw_buffer,
                        std::complex<float>* visibility_buffer,
                        const aocommon::BandData& cur_band,
                        float* weight_buffer, std::complex<float>* model_buffer,
                        const bool* is_selected);

  // Equivalent functionality to 'ApplyWeightsAndCorrections' however written so
  // that it can apply the corrections inlines and out of order e.g. from within
  // a gridder that wants to get them "on the fly" As opposed to sequentially
  // being called in order one row at a time like 'ApplyWeightsAndCorrections'
  // is
  template <size_t PolarizationCount, GainMode GainEntry>
  inline std::complex<float>& InlineApplyWeightsAndCorrections(
      const std::vector<std::string>& antenna_names,
      const aocommon::BandData& curBand, float* scratch_weight_buffer,
      float* weight_buffer, std::complex<float>* visibility_buffer,
      size_t visibility_index, size_t row, size_t channel,
      const MSProvider::MetaData& metaData) {
    const size_t channel_count = curBand.ChannelCount();
    const size_t data_size = PolarizationCount * channel_count;

    // Calculate a row at a time, and cache up to 256 rows before throwing one
    // away (well overwriting instead) This avoids excessive duplicate
    // computation.
    thread_local static std::vector<aocommon::UVector<std::complex<float>>>
        visibility_row(256, aocommon::UVector<std::complex<float>>(data_size));
    thread_local static size_t last_row = std::numeric_limits<size_t>::max();
    thread_local static uint8_t row_index = 0;

    // We might be applying these same weights/corrections multiple times, so
    // ensure it only counts towards the totals the first time around
    bool already_counted = already_counted_weights[row];
    if (!already_counted) {
      already_counted_weights[row] = true;
    }

    if (last_row != row) {
      last_row = row;
      ++row_index;

      // Read the visibilities into our cache so we can modify them in place
      // without touching the original underlying data
      std::copy_n(&visibility_buffer[row * data_size], data_size,
                  visibility_row[row_index].data());
      float* row_scratch_weights = &scratch_weight_buffer[row * channel_count];
      float* row_weight_iter = &weight_buffer[row * data_size];
      std::complex<float>* row_visibility_iter =
          visibility_row[row_index].data();

      // Apply corrections
      if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
        const bool apply_forward = GetPsfMode() == PsfMode::kDirectionDependent;
        const bool apply_beam =
            settings_.applyFacetBeam || settings_.gridWithBeam;
        if (apply_beam && visibility_modifier_.HasH5Parm()) {
          // This has been temporarily removed for prototyping purposes. We will
          // need to add it back once we have proven the prototype.
          exit(721);
        } else if (apply_beam) {
          // This has been temporarily removed for prototyping purposes. We will
          // need to add it back once we have proven the prototype.
          exit(722);
        } else if (visibility_modifier_.HasH5Parm()) {
          visibility_modifier_
              .ApplyConjugatedParmResponse<PolarizationCount, GainEntry>(
                  visibility_row[row_index].data(), row_weight_iter,
                  row_scratch_weights, ms_index_, channel_count,
                  antenna_names.size(), metaData.antenna1, metaData.antenna2,
                  apply_forward, time_offsets_[ms_index_][row],
                  !already_counted);
        }
      }

      // Adjust these values locally then only set them once to the shared
      // variable, this way we reduce contention and reduce the need for locking
      double local_visibility_weight_sum = 0.0;
      double local_max_gridded_weight = 0.0;
      double local_total_weight = 0.0;
      size_t local_visibility_count = 0;
      // Apply visibility and imaging weights
      for (size_t ch = 0; ch != channel_count; ++ch) {
        for (size_t p = 0; p != PolarizationCount; ++p) {
          double cum_weight = *row_weight_iter * row_scratch_weights[ch];

          if (p == 0 && cum_weight != 0.0) {
            // Visibility weight sum is the sum of weights excluding imaging
            // weights
            local_visibility_weight_sum += *row_weight_iter;
            local_max_gridded_weight =
                std::max(cum_weight, local_max_gridded_weight);
            ++local_visibility_count;
          }
          // Total weight includes imaging weights
          local_total_weight += cum_weight;

          *row_visibility_iter *= cum_weight;
          ++row_visibility_iter;
          ++row_weight_iter;
        }
      }
      total_weight_mutex_.lock();
      if (!already_counted) {
        visibility_weight_sum_ += local_visibility_weight_sum;
        max_gridded_weight_ =
            std::max(local_max_gridded_weight, max_gridded_weight_);
        gridded_visibility_count_ += local_visibility_count;
        total_weight_ += local_total_weight;
      }
      total_weight_mutex_.unlock();
    }

    // Return the actual visibility that has been asked for from the cache
    return visibility_row[row_index][channel];
  }

  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyWeightsAndCorrections(const std::vector<std::string>& antenna_names,
                                  InversionRow& row_data,
                                  const aocommon::BandData& cur_band,
                                  float* weight_buffer,
                                  const MSProvider::MetaData& meta_data);

  template <size_t PolarizationCount, GainMode GainEntry>
  void WriteInstrumentalVisibilities(
      MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
      const aocommon::BandData& cur_band, std::complex<float>* buffer,
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
  /// MeasurementSetCount() determine the index in the _writerGroupLocks vector,
  /// having size FacetGroupCount() * MeasurementSetCount(). These variable are
  /// only relevant for prediction.
  size_t facet_group_index_ = 0;
  size_t ms_index_ = 0;
  /// @see SetIsFacet()
  bool is_facet_ = false;
  double image_padding_ = 1.0;
  size_t image_width_ = 0.0;
  size_t image_height_ = 0.0;
  size_t trim_width_ = 0;
  size_t trim_height_ = 0;
  size_t w_grid_size_ = 0;
  size_t actual_w_grid_size_ = 0;
  std::vector<std::unique_ptr<MSProvider>> measurement_sets_;
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
  GainMode gain_mode_ = GainMode::kDiagonal;
  bool is_complex_ = false;
  WeightMode weighting_ = WeightMode(WeightMode::UniformWeighted);
  bool is_first_task_ = false;
  std::vector<MSSelection> selections_;
  VisibilityWeightingMode visibility_weighting_mode_ =
      VisibilityWeightingMode::NormalVisibilityWeighting;
  GriddingKernelMode grid_mode_ = GriddingKernelMode::KaiserBessel;
  bool store_imaging_weights_ = false;
  double theoretical_beam_size_ = 0.0;

  bool has_frequencies_ = false;
  double freq_high_ = 0.0;
  double freq_low_ = 0.0;
  double band_start_ = 0.0;
  double band_end_ = 0.0;
  double start_time_ = 0.0;

 public:
  std::vector<bool> already_counted_weights;
  std::map<size_t, std::vector<size_t>> time_offsets_;

  size_t gridded_visibility_count_ = 0;
  std::mutex total_weight_mutex_;
  double total_weight_ = 0.0;
  double max_gridded_weight_ = 0.0;
  double visibility_weight_sum_ = 0.0;

  aocommon::UVector<float> scratch_image_weights_;

 private:
  /// Initialized in StartMeasurementSet(), used in WriteCollapsedVisibilities()
  /// to expand visibilities into.
  aocommon::UVector<std::complex<float>> scratch_model_data_;

  std::unique_ptr<MSReader> predict_reader_;
  GriddingTaskManager* writer_lock_manager_ = nullptr;

  VisibilityModifier visibility_modifier_;
};

template <size_t PolarizationCount>
inline std::complex<float>& MSGridderBase::InlineApplyWeightsAndCorrections(
    const std::vector<std::string>& antenna_names,
    const aocommon::BandData& cur_band, float* scratch_weight_buffer,
    float* weight_buffer, std::complex<float>* visibility_buffer,
    size_t visibility_index, size_t row, size_t channel,
    const MSProvider::MetaData& meta_data) {
  switch (gain_mode_) {
    case GainMode::kXX:
      if constexpr (PolarizationCount == 1) {
        return InlineApplyWeightsAndCorrections<PolarizationCount,
                                                GainMode::kXX>(
            antenna_names, cur_band, scratch_weight_buffer, weight_buffer,
            visibility_buffer, visibility_index, row, channel, meta_data);
      }
      break;
    case GainMode::kYY:
      if constexpr (PolarizationCount == 1) {
        return InlineApplyWeightsAndCorrections<PolarizationCount,
                                                GainMode::kYY>(
            antenna_names, cur_band, scratch_weight_buffer, weight_buffer,
            visibility_buffer, visibility_index, row, channel, meta_data);
      }
      break;
    case GainMode::kDiagonal:
      if constexpr (PolarizationCount == 1 || PolarizationCount == 2) {
        return InlineApplyWeightsAndCorrections<PolarizationCount,
                                                GainMode::kDiagonal>(
            antenna_names, cur_band, scratch_weight_buffer, weight_buffer,
            visibility_buffer, visibility_index, row, channel, meta_data);
      }
      break;
    case GainMode::kFull:
      if constexpr (PolarizationCount == 4) {
        return InlineApplyWeightsAndCorrections<PolarizationCount,
                                                GainMode::kFull>(
            antenna_names, cur_band, scratch_weight_buffer, weight_buffer,
            visibility_buffer, visibility_index, row, channel, meta_data);
      } else {
        throw std::runtime_error(
            "Invalid combination of visibility polarizations and gain mode");
      }
      break;
    default:
      throw std::runtime_error("Unknown gain mode");
  }
}

template <size_t PolarizationCount>
inline void MSGridderBase::ApplyWeightsAndCorrections(
    const std::vector<std::string>& antenna_names, InversionRow& row_data,
    const aocommon::BandData& cur_band, float* weight_buffer,
    const MSProvider::MetaData& meta_data) {
  switch (gain_mode_) {
    case GainMode::kXX:
      if constexpr (PolarizationCount == 1) {
        ApplyWeightsAndCorrections<PolarizationCount, GainMode::kXX>(
            antenna_names, row_data, cur_band, weight_buffer, meta_data);
      }
      break;
    case GainMode::kYY:
      if constexpr (PolarizationCount == 1) {
        ApplyWeightsAndCorrections<PolarizationCount, GainMode::kYY>(
            antenna_names, row_data, cur_band, weight_buffer, meta_data);
      }
      break;
    case GainMode::kDiagonal:
      if constexpr (PolarizationCount == 1 || PolarizationCount == 2) {
        ApplyWeightsAndCorrections<PolarizationCount, GainMode::kDiagonal>(
            antenna_names, row_data, cur_band, weight_buffer, meta_data);
      }
      break;
    case GainMode::kFull:
      if constexpr (PolarizationCount == 4) {
        ApplyWeightsAndCorrections<PolarizationCount, GainMode::kFull>(
            antenna_names, row_data, cur_band, weight_buffer, meta_data);
      } else {
        throw std::runtime_error(
            "Invalid combination of visibility polarizations and gain mode");
      }
      break;
  }
}

template <size_t PolarizationCount>
inline void MSGridderBase::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& cur_band, std::complex<float>* buffer,
    MSProvider::MetaData& meta_data) {
  switch (gain_mode_) {
    case GainMode::kXX:
      if constexpr (PolarizationCount == 1) {
        WriteInstrumentalVisibilities<PolarizationCount, GainMode::kXX>(
            ms_provider, antenna_names, cur_band, buffer, meta_data);
      }
      break;
    case GainMode::kYY:
      if constexpr (PolarizationCount == 1) {
        WriteInstrumentalVisibilities<PolarizationCount, GainMode::kYY>(
            ms_provider, antenna_names, cur_band, buffer, meta_data);
      }
      break;
    case GainMode::kDiagonal:
      if constexpr (PolarizationCount == 1 || PolarizationCount == 2) {
        WriteInstrumentalVisibilities<PolarizationCount, GainMode::kDiagonal>(
            ms_provider, antenna_names, cur_band, buffer, meta_data);
      }
      break;
    case GainMode::kFull:
      if constexpr (PolarizationCount == 4)
        WriteInstrumentalVisibilities<PolarizationCount, GainMode::kFull>(
            ms_provider, antenna_names, cur_band, buffer, meta_data);
      else
        throw std::runtime_error(
            "Invalid combination of visibility polarizations and gain mode");
      break;
  }
}

#endif
