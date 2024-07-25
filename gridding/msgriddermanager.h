#ifndef WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
#define WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_

#include <mutex>
#include <vector>

#include "h5solutiondata.h"
#include "msgridderbase.h"

#include "../main/settings.h"
#include "../scheduling/griddingresult.h"
#include "../scheduling/griddingtask.h"
#include "../structures/resources.h"

namespace wsclean {

class GriddingTaskManager;

/**
 * The MSGridderManager is a middle layer between GriddingTaskManager and
 * MSGridderBase derived classes.
 *
 * GriddingTaskManager is solely responsible for scheduling MSGridderBase
 * derived classes are responsible for gridding (inversion/predict)
 *
 * MSGridderManager is responsible for:
 *   1. Managing gridding specific data thats lifecycle/scope is larger than a
 *      single gridder.
 *      e.g. h5parm solution data. Reading a single visibility buffer that
 *      multiple gridders will all operate on.
 *   2. Setting up requirements of the gridder classes.
 *      e.g. initialising MS readers
 *   3. Initialising gridder classes
 *   4. Starting gridder classes
 *   5. Performing final processing of the results before passing them back
 *      through the task manager
 */
class MSGridderManager {
 public:
  MSGridderManager(const Settings& settings,
                   const H5SolutionData& solution_data)
      : settings_(settings), solution_data_(solution_data) {}
  ~MSGridderManager() {}

  MSGridderManager(const MSGridderManager&) = delete;
  MSGridderManager& operator=(const MSGridderManager&) = delete;

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

  /** Initializes 'gridder' with values that are equal for all facets. */
  void InitializeGridderForTask(MSGridderBase& gridder,
                                const GriddingTask& task,
                                GriddingTaskManager* writer_lock_manager);

  /** Initializes 'gridder' with facet-specific values. */
  void InitializeGridderForFacet(MSGridderBase& gridder,
                                 GriddingTask::FacetData& facet_task);

  const Settings& settings_;
  const H5SolutionData& solution_data_;
};

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
