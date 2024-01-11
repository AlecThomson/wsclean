#ifndef WSCLEAN_MSHELPER_H_
#define WSCLEAN_MSHELPER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>

#include "../msproviders/msdatadescription.h"
#include "../msproviders/partitionedms.h"
#include "../structures/imagingtableentry.h"
#include "../structures/msselection.h"

#include "settings.h"

/**
 * Class with helper routines for managing measurement sets.
 */
class MsHelper {
 public:
  explicit MsHelper(
      const Settings& settings, const MSSelection& global_selection,
      const std::vector<aocommon::MultiBandData>& ms_bands,
      const std::vector<PartitionedMS::Handle>& partitioned_ms_handles)
      : settings_(settings),
        global_selection_(global_selection),
        ms_bands_(ms_bands),
        partitioned_ms_handles_(partitioned_ms_handles) {}

  std::vector<std::unique_ptr<MSDataDescription>> InitializeMsList(
      const ImagingTableEntry& entry) const;

 private:
  const Settings& settings_;
  const MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  const std::vector<PartitionedMS::Handle>& partitioned_ms_handles_;
};

#endif