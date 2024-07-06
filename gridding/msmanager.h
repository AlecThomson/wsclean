#ifndef MS_MANAGER_H
#define MS_MANAGER_H

#include <mutex>
#include <vector>

#include "h5solutiondata.h"
#include "msgridderbase.h"

#include "../main/settings.h"
#include "../scheduling/griddingresult.h"
#include "../scheduling/griddingtask.h"
#include "../structures/resources.h"

class MSGridderBase;

class MSManager {
 public:
  struct Data {
   public:
    MSProvider* ms_provider = nullptr;
    size_t internal_ms_index = 0;
    size_t original_ms_index = 0;
    size_t dataDescId = 0;
    aocommon::BandData bandData;
    size_t startChannel = 0;
    size_t endChannel = 0;
    mutable size_t matchingRows = 0;
    mutable size_t totalRowsProcessed = 0;
    size_t rowStart = 0;
    size_t rowEnd = 0;
    std::vector<std::string> antenna_names;
    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(bandData, startChannel, endChannel);
    }
    void InitializeBandData(const casacore::MeasurementSet& ms,
                            const MSSelection& selection);
  };
  struct FacetData {
   public:
    double minW = 0.0;
    double maxW = 0.0;
    double maxWWithFlags = 0.0;
    double maxBaselineUVW = 0.0;
    double maxBaselineInM = 0.0;
    double integrationTime = 0.0;
  };

  MSProvider& MeasurementSet(size_t internal_ms_index) const {
    return *std::get<0>(measurement_sets_[internal_ms_index]);
  }
  const MSSelection& Selection(size_t internal_ms_index) const {
    return std::get<1>(measurement_sets_[internal_ms_index]);
  }
  /**
   * Maps an internal ms index to its original (command line) index.
   */
  size_t Index(size_t internal_ms_index) const {
    return std::get<2>(measurement_sets_[internal_ms_index]);
  }
  size_t Count() const { return measurement_sets_.size(); }

  /**
   * @param index the original command line index of this measurement set. This
   * allows associating the measurement set with the right h5parm solution file.
   */
  void Add(std::unique_ptr<MSProvider> ms_provider,
           const MSSelection& selection, size_t index) {
    measurement_sets_.emplace_back(std::move(ms_provider), selection, index);
  }

  double StartTime() const { return ms_limits_.start_time_; }

  void InitializeMSDataVector(const std::vector<MSGridderBase*>& gridders);

  // Computes/stores per MS metadata that is relevant/shared for all
  // facets/gridders
  std::vector<Data> ms_data_vector_;
  // Computes/stores per MS metadata that is unique to an
  // individrelevant/shared for all facets/gridders
  /// First index is gridder/facet, second index is MS
  std::vector<std::vector<FacetData>> ms_facet_data_vector_;

 private:
  template <size_t NPolInMSProvider>
  void CalculateWLimits(FacetData& ms_facet_data, Data& ms_data,
                        MSGridderBase* gridder);

  static std::vector<std::string> GetAntennaNames(
      const casacore::MSAntenna& antenna);

  void InitializeMeasurementSet(
      Data& ms_data, std::vector<std::vector<FacetData>>& ms_facet_data_vector,
      const std::vector<MSGridderBase*>& gridders, bool is_cached);

  // Computes/stores frequency limits across all measurement sets
  struct {
   public:
    bool has_frequencies_ = false;
    double freq_high_ = 0.0;
    double freq_low_ = 0.0;
    double band_start_ = 0.0;
    double band_end_ = 0.0;
    double start_time_ = 0.0;

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
  } ms_limits_;

  /// tuple consists of the ms, selection, and index.
  std::vector<std::tuple<std::unique_ptr<MSProvider>, MSSelection, size_t>>
      measurement_sets_;
};

#endif
