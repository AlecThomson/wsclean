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
    size_t data_desc_id = 0;
    aocommon::BandData band_data;
    size_t start_channel = 0;
    size_t end_channel = 0;
    size_t matching_rows = 0;
    size_t total_rows_processed = 0;
    size_t row_start = 0;
    size_t row_end = 0;

    double min_w = 0.0;
    double max_w = 0.0;
    double max_w_with_flags = 0.0;
    double max_baseline_uvw = 0.0;
    double max_baseline_meters = 0.0;
    double integration_time = 0.0;

    std::vector<std::string> antenna_names;
    std::shared_ptr<std::vector<double>> unique_times;

    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(band_data, start_channel, end_channel);
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

  double StartTime() const { return ms_limits_.start_time; }

  void InitializeMS();
  void InitializeMSDataVector(const std::vector<MSGridderBase*>& gridders,
                              double w_limit, bool has_solution_data);

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
                                bool is_cached, bool has_solution_data);

  /** Computes/stores limits across all measurement sets */
  struct {
   public:
    bool has_frequencies = false;
    double highest_frequency = 0.0;
    double lowest_frequency = 0.0;
    double band_start = 0.0;
    double band_end = 0.0;
    double start_time = 0.0;

    double max_w = 0.0;
    double min_w = 0.0;
    double max_baseline = 0.0;

    void Calculate(const aocommon::BandData& selected_band,
                   double _start_time) {
      if (has_frequencies) {
        lowest_frequency =
            std::min(lowest_frequency, selected_band.LowestFrequency());
        highest_frequency =
            std::max(highest_frequency, selected_band.HighestFrequency());
        band_start = std::min(band_start, selected_band.BandStart());
        band_end = std::max(band_end, selected_band.BandEnd());
        start_time = std::min(start_time, _start_time);
      } else {
        lowest_frequency = selected_band.LowestFrequency();
        highest_frequency = selected_band.HighestFrequency();
        band_start = selected_band.BandStart();
        band_end = selected_band.BandEnd();
        start_time = _start_time;
        has_frequencies = true;
      }
    }
    void Validate() {
      if (min_w > max_w) {
        min_w = max_w;
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
