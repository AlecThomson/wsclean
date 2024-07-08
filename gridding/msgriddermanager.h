#ifndef MS_GRIDDER_MANAGER_H
#define MS_GRIDDER_MANAGER_H

#include <mutex>
#include <vector>

#include "h5solutiondata.h"
#include "msgridderbase.h"
#include "msmanager.h"

#include "../main/settings.h"
#include "../scheduling/griddingresult.h"
#include "../scheduling/griddingtask.h"
#include "../structures/resources.h"

class GriddingTaskManager;

/**
 * The MSGridderManager is a middle layer between GriddingTaskManager and
 * MSGridderBase derived classes.
 *
 * GriddingTaskManager is solely responsible for scheduling
 * MSGridderBase derived classes are responsible for gridding
 * (inversion/predict) MSGridderManager is responsible for: Managing gridding
 * specific data thats lifecycle/scope is larger than a single gridder. e.g.
 * h5parm solution data or reading a single visibility buffer that multiple
 * gridders will all operate on Setting up requirements of the gridder classes
 * e.g. initialising MS readers Initialising gridder classes Starting gridder
 * classes Performing final processing of the results before passing them back
 * through the task manager
 */
class MSGridderManager {
 public:
  MSGridderManager(const Settings& settings,
                   const H5SolutionData& solution_data)
      : settings_(settings), solution_data_(solution_data) {}
  ~MSGridderManager() {}

  MSGridderManager(const MSGridderManager&) = delete;
  MSGridderManager& operator=(const MSGridderManager&) = delete;

  void InitializeMS(GriddingTask& task);

  void InitializeGridders(GriddingTask& task,
                          const std::vector<size_t>& facet_indices,
                          const Resources& resources,
                          std::vector<GriddingResult::FacetData>& facet_results,
                          GriddingTaskManager* writer_lock_manager);
  void Invert();
  void Predict();
  void ProcessResults(std::mutex& result_mutex, GriddingResult& result,
                      bool store_common_info);

 private:
  std::unique_ptr<MSGridderBase> ConstructGridder(
      const Resources& resources) const;
  struct GriddingFacetTask {
    std::unique_ptr<MSGridderBase> facet_gridder;
    GriddingTask::FacetData& facet_task;
    GriddingResult::FacetData& facet_result;
  };
  std::vector<GriddingFacetTask> facet_tasks_;

  inline void InitializeMSDataVectors() {
    std::vector<MSGridderBase*> gridders;
    gridders.reserve(facet_tasks_.size());
    for (auto& [gridder, facet_task, facet_result] : facet_tasks_) {
      gridders.push_back(gridder.get());
    }
    measurement_sets_.InitializeMSDataVector(gridders);
  }

  /** Initializes 'gridder' with values that are equal for all facets. */
  void InitializeGridderForTask(MSGridderBase& gridder,
                                const GriddingTask& task,
                                GriddingTaskManager* writer_lock_manager);

  /** Initializes 'gridder' with facet-specific values. */
  void InitializeGridderForFacet(MSGridderBase& gridder,
                                 GriddingTask::FacetData& facet_task);

  const Settings& settings_;
  const H5SolutionData& solution_data_;
  MSManager measurement_sets_;
};
#endif
