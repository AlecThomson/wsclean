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

  /**
   * Creates gridding / inversion tasks for PSF images.
   * @param facet_group Group with all imaging table entries (facets) for a
   * single PSF image.
   * When faceting is disabled, the group contains only one entry.
   * @param combine_facets True: Create a single task for all facets.
   * False: Create an individual task for each facet.
   */
  std::vector<GriddingTask> CreatePsfTasks(
      const ImagingTable::Group& facet_group,
      ImageWeightCache& image_weight_cache, bool combine_facets,
      bool is_first_task);

  /**
   * Creates gridding / inversion tasks.
   * @param facet_group Group with all imaging table entries (facets) for a
   * single PSF image.
   * When faceting is disabled, the group contains only one entry.
   * @param is_first_task True for the first gridding task. WSClean produces
   * more verbose logging output for the first task.
   * @param is_first_inversion True if the task is part of the first set of
   * gridding / inversion tasks.
   * @param combine_facets True: Create a single task for all facets.
   * False: Create an individual task for each facet.
   */
  std::vector<GriddingTask> CreateInvertTasks(
      const ImagingTable::Group& facet_group,
      ImageWeightCache& image_weight_cache, bool combine_facets,
      bool is_first_task, bool is_first_inversion,
      std::vector<std::unique_ptr<AverageBeam>>&& average_beams);

  /**
   * Creates a degridding / predict task.
   */
  GriddingTask CreatePredictTask(const ImagingTableEntry& entry,
                                 ImageWeightCache& image_weight_cache,
                                 std::vector<aocommon::Image>&& model_images,
                                 std::unique_ptr<AverageBeam> average_beam);

  const std::vector<std::unique_ptr<MetaDataCache>>& GetMetaDataCache() const {
    return meta_data_cache_;
  }

  void SetMetaDataCacheEntry(const ImagingTableEntry& entry,
                             std::unique_ptr<MetaDataCache> cache) {
    assert(entry.index < meta_data_cache_.size());
    meta_data_cache_[entry.index] = std::move(cache);
  }

  long double GetFacetCorrectionFactor(const ImagingTableEntry& entry) const {
    return meta_data_cache_[entry.index]->correctionSum / entry.imageWeight;
  }

 private:
  /// Add a facet structure to an existing task.
  /// @param task An existing task, which may not have FacetData yet.
  /// @param entry ImagingTableEntry with facet-specific data.
  /// @param average_beam Average beam value for the facet. May be empty.
  void AddFacet(GriddingTask& task, const ImagingTableEntry& entry,
                std::unique_ptr<AverageBeam> average_beam = nullptr);

  /// Determine the polarization value for an invert or predict task.
  aocommon::PolarizationEnum DeterminePolarization(
      const ImagingTableEntry& entry) const;

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