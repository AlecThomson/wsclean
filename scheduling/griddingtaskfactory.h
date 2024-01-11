#ifndef GRIDDING_TASK_FACTORY_H_
#define GRIDDING_TASK_FACTORY_H_

#include <memory>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include "../main/imageweightinitializer.h"
#include "../main/mshelper.h"
#include "../structures/imagingtableentry.h"
#include "../structures/observationinfo.h"

#include "griddingtask.h"
#include "metadatacache.h"

class GriddingTaskFactory {
 public:
  explicit GriddingTaskFactory(const MsHelper& ms_helper,
                               const ImageWeightInitializer& initializer,
                               const ObservationInfo& observation_info,
                               double l_shift, double m_shift,
                               size_t imaging_table_size)
      : ms_helper_(ms_helper),
        image_weight_initializer_{initializer},
        observation_info_{observation_info},
        l_shift_{l_shift},
        m_shift_{m_shift},
        meta_data_cache_{imaging_table_size} {}

  GriddingTask CreatePsfTask(const ImagingTableEntry& entry,
                             ImageWeightCache& image_weight_cache,
                             bool is_first_task);

  /**
   * Creates a degridding / inversion task.
   * @param is_first_task True for the first gridding task. WSClean produces
   * more verbose logging output for the first task.
   * @param is_first_inversion True if the task is part of the first set of
   * degridding / inversion tasks.
   */
  GriddingTask CreateInvertTask(const ImagingTableEntry& entry,
                                ImageWeightCache& image_weight_cache,
                                bool is_first_task, bool is_first_inversion,
                                std::unique_ptr<AverageBeam> average_beam);

  GriddingTask CreatePredictTask(const ImagingTableEntry& entry,
                                 ImageWeightCache& image_weight_cache,
                                 std::vector<aocommon::Image>&& model_images,
                                 std::unique_ptr<AverageBeam> average_beam);

  const std::vector<std::unique_ptr<MetaDataCache>>& GetMetaDataCache() const {
    return meta_data_cache_;
  }

  void SetMetaDataCacheEntry(const ImagingTableEntry& entry,
                             std::unique_ptr<MetaDataCache> cache) {
    meta_data_cache_[entry.index] = std::move(cache);
  }

  long double GetFacetCorrectionFactor(const ImagingTableEntry& entry) const {
    const std::unique_ptr<MetaDataCache>& cache_entry =
        meta_data_cache_[entry.index];
    return cache_entry->correctionSum / entry.imageWeight;
  }

 private:
  /// Common code for initializing a task.
  /// After this function, the remaining task fields should be initialized.
  GriddingTask CreateBase(const ImagingTableEntry& entry,
                          ImageWeightCache& image_weight_cache,
                          bool is_first_task);

 private:
  const MsHelper& ms_helper_;
  const ImageWeightInitializer& image_weight_initializer_;
  const ObservationInfo& observation_info_;
  const double l_shift_, m_shift_;
  std::vector<std::unique_ptr<MetaDataCache>> meta_data_cache_;
};

#endif