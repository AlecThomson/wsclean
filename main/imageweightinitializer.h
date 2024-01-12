#ifndef WSCLEAN_INITIALIZER_H_
#define WSCLEAN_INITIALIZER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include "../io/imageweightcache.h"
#include "../msproviders/msdatadescription.h"
#include "../msproviders/partitionedms.h"
#include "../structures/imageweights.h"
#include "../structures/imagingtable.h"
#include "../structures/msselection.h"

#include "settings.h"

/**
 * Helper class with various initialization routines for WSClean.
 */
class ImageWeightInitializer {
 public:
  explicit ImageWeightInitializer(
      const Settings& settings, const MSSelection& global_selection,
      const std::vector<aocommon::MultiBandData>& ms_bands,
      const std::vector<PartitionedMS::Handle>& reordered_ms_handles)
      : settings_(settings),
        global_selection_(global_selection),
        ms_bands_(ms_bands),
        reordered_ms_handles_(reordered_ms_handles) {}

  const Settings& GetSettings() const { return settings_; }

  std::shared_ptr<ImageWeights> Initialize(
      const ImagingTableEntry& entry,
      const std::vector<std::unique_ptr<MSDataDescription>>& ms_list,
      ImageWeightCache& cache) const;

  void InitializeMf(const ImagingTable& imaging_table, ImageWeightCache& cache);

 private:
  const Settings& settings_;
  const MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  const std::vector<PartitionedMS::Handle>& reordered_ms_handles_;
};

#endif