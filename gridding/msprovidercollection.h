#ifndef WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_
#define WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_

#include <vector>

#include <aocommon/logger.h>

#include "../structures/imageweights.h"
#include "../msproviders/msprovider.h"
#include "../structures/msselection.h"

namespace wsclean {

class MSGridderBase;

/**
 * MsProviderCollection is responsible for managing the set of measurement sets
 * on which one or multiple gridders will operate, as well as computing and
 * storing metadata about the MS that will be used by the gridders
 *
 * Metadata is split and stored in two different structs:
 *  1. `MsData` contains per MS data that is shared across all facets in a
 *     facet group. One instance per MS
 *  2. `Limits` contains data that is calculated across all MSW. One instance
 *     for multiple MS
 */
class MsProviderCollection {
 public:
  struct MsData {
   public:
    MSProvider* ms_provider = nullptr;
    size_t internal_ms_index = 0;
    size_t original_ms_index = 0;
    size_t dataDescId = 0;
    aocommon::BandData bandData;
    size_t startChannel = 0;
    size_t endChannel = 0;
    size_t matchingRows = 0;
    // TODO: remove mutable once we have restructured data reading out of the
    // gridders and into the gridder manager
    mutable size_t totalRowsProcessed = 0;
    size_t rowStart = 0;
    size_t rowEnd = 0;

    double minW = 0.0;
    double maxW = 0.0;
    double maxWWithFlags = 0.0;
    double maxBaselineUVW = 0.0;
    double maxBaselineInM = 0.0;
    double integrationTime = 0.0;

    std::vector<std::string> antenna_names;
    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(bandData, startChannel, endChannel);
    }
    void InitializeBandData(const casacore::MeasurementSet& ms,
                            const MSSelection& selection);
  };

  MSProvider& MeasurementSet(size_t internal_ms_index) const {
    return *std::get<0>(ms_provider_collection_[internal_ms_index]);
  }
  const MSSelection& Selection(size_t internal_ms_index) const {
    return std::get<1>(ms_provider_collection_[internal_ms_index]);
  }
  /**
   * Maps an internal ms index to its original (command line) index.
   */
  size_t Index(size_t internal_ms_index) const {
    return std::get<2>(ms_provider_collection_[internal_ms_index]);
  }
  size_t Count() const { return ms_provider_collection_.size(); }

  /**
   * @param index the original command line index of this measurement set. This
   * allows associating the measurement set with the right h5parm solution file.
   */
  void Add(std::unique_ptr<MSProvider> ms_provider,
           const MSSelection& selection, size_t index) {
    ms_provider_collection_.emplace_back(std::move(ms_provider), selection,
                                         index);
  }

  double StartTime() const { return ms_limits_.start_time_; }

  void InitializeMS();
  void InitializeMSDataVector(const std::vector<MSGridderBase*>& gridders,
                              double w_limit);

  /** Computes/stores per MS metadata that is shared for all facets/gridders in
   * a facet group */
  std::vector<MsData> ms_data_vector_;

 private:
  template <size_t NPolInMSProvider>
  void CalculateMsLimits(MsData& ms_data, double pixel_size_x,
                         double pixel_size_y, size_t image_width,
                         size_t image_height,
                         const ImageWeights* image_weights);

  static std::vector<std::string> GetAntennaNames(
      const casacore::MSAntenna& antenna);

  void InitializeMeasurementSet(MsData& ms_data,
                                const std::vector<MSGridderBase*>& gridders,
                                bool is_cached);

  /** Computes/stores limits across all measurement sets */
  struct {
   public:
    bool has_frequencies_ = false;
    double freq_high_ = 0.0;
    double freq_low_ = 0.0;
    double band_start_ = 0.0;
    double band_end_ = 0.0;
    double start_time_ = 0.0;

    double max_w_ = 0.0;
    double min_w_ = 0.0;
    double maxBaseline_ = 0.0;

    void Calculate(const aocommon::BandData& selectedBand, double startTime) {
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
    void Validate() {
      if (min_w_ > max_w_) {
        min_w_ = max_w_;
        aocommon::Logger::Error
            << "*** Error! ***\n"
               "*** Calculating maximum and minimum w values failed! Make sure "
               "the "
               "data selection and scale settings are correct!\n"
               "***\n";
      }
    }
  } ms_limits_;

  /// tuple consists of the ms, selection, and index.
  std::vector<std::tuple<std::unique_ptr<MSProvider>, MSSelection, size_t>>
      ms_provider_collection_;
};

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_
