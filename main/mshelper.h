#ifndef WSCLEAN_MSHELPER_H_
#define WSCLEAN_MSHELPER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>

#include "../msproviders/msdatadescription.h"
#include "../msproviders/reorderedms.h"
#include "../structures/imagingtable.h"
#include "../structures/msselection.h"

#include "settings.h"

/**
 * Class with helper routines for managing measurement sets.
 */
class MsHelper {
 public:
  explicit MsHelper(const Settings& settings,
                    const MSSelection& global_selection,
                    const std::vector<aocommon::MultiBandData>& ms_bands)
      : settings_{settings},
        global_selection_{global_selection},
        ms_bands_{ms_bands},
        reordered_ms_handles_{} {}

  const std::vector<ReorderedMs::Handle>& GetReorderedMsHandles() const {
    return reordered_ms_handles_;
  }

  void PerformReordering(const ImagingTable& imaging_table,
                         bool is_predict_mode);

  std::vector<std::unique_ptr<MSDataDescription>> InitializeMsList(
      const ImagingTableEntry& entry) const;

 private:
  const Settings& settings_;
  const MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  std::vector<ReorderedMs::Handle> reordered_ms_handles_;
};

#endif