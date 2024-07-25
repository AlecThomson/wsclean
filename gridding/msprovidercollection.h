#ifndef WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_
#define WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_

#include <vector>

#include "../msproviders/msprovider.h"
#include "../structures/msselection.h"

namespace wsclean {

class MSGridderBase;

/**
 * MsProviderCollection is responsible for managing the set of measurement sets
 * on which one or multiple gridders will operate, as well as computing and
 * storing metadata about the MS that will be used by the gridders
 *
 * Metadata is split and stored in three different structs:
 *  1. `MsData` contains shared data for a MS that is the same across all facets
 *     in a facet group (One instance per MS)
 *  2. `FacetData` contains facet specific data for a MS that is unique to each
 *     facet and can't be shared across other facets in the facet group (Instace
 *     per facet per MS)
 *  3. `Limits` contains data that is calculated across all MS and is (One
 *     instance for multiple MS)
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

  void InitializeMS(size_t num_facets);
  void InitializeMSDataVector(const std::vector<MSGridderBase*>& gridders);

  /** Computes/stores per MS metadata that is shared for all facets/gridders in
   * a facet group */
  std::vector<MsData> ms_data_vector_;

  /** Computes and stores per MS metadata that is unique to each individual
   * facet/gridder in a facet group.
   * Outer index is facet/gridder.
   * Inner index is MS.
   */
  std::vector<std::vector<FacetData>> ms_facet_data_vector_;

 private:
  template <size_t NPolInMSProvider>
  void CalculateWLimits(FacetData& ms_facet_data, MsData& ms_data,
                        MSGridderBase* gridder);

  static std::vector<std::string> GetAntennaNames(
      const casacore::MSAntenna& antenna);

  void InitializeMeasurementSet(
      MsData& ms_data,
      std::vector<std::vector<FacetData>>& ms_facet_data_vector,
      const std::vector<MSGridderBase*>& gridders, bool is_cached);

  /** Computes/stores frequency limits across all measurement sets */
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
      ms_provider_collection_;
};

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_PROVIDER_COLLECTION_H_
