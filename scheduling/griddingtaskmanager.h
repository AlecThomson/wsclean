#ifndef GRIDDING_TASK_MANAGER_H_
#define GRIDDING_TASK_MANAGER_H_

#include <cstring>
#include <functional>
#include <memory>
#include <vector>

#include <aocommon/lane.h>
#include <aocommon/polarization.h>

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"
#include "../structures/msselection.h"

#include "../msproviders/msdatadescription.h"

#include "griddingtask.h"
#include "griddingresult.h"

class MSGridderBase;
class Resources;
class Settings;

class GriddingTaskManager {
 public:
  class WriterLock {
   public:
    virtual ~WriterLock() = default;
  };

 public:
  explicit GriddingTaskManager(const Settings& settings);

  virtual ~GriddingTaskManager();

  const Settings& GetSettings() const { return settings_; }

  void SetWriterLockManager(GriddingTaskManager& manager) {
    writer_lock_manager_ = &manager;
  }

  /**
   * Initialize writer groups. Call this function before scheduling Predict
   * tasks in order to initialize the writer locks.
   *
   * @param nWriterGroups The number of writer groups.
   */
  virtual void Start([[maybe_unused]] size_t nWriterGroups) {}

  /**
   * Obtains a lock for the given @p writer_group_index.
   * The default implementation returns a dummy lock: Since it runs all tasks
   * sequentially, locking is not needed.
   * @return A lock object. Destroying the object releases the lock.
   */
  virtual std::unique_ptr<WriterLock> GetLock(
      [[maybe_unused]] size_t writer_group_index) {
    return nullptr;
  }

  /**
   * Add the given task to the queue of tasks to be run. After finishing
   * the task, the callback is called with the results. The callback will
   * always run in the thread of the caller.
   * Depending on the type of gridding task manager, this call might block.
   *
   * This implementation runs the task directly and blocks until done.
   */
  virtual void Run(GriddingTask&& task,
                   std::function<void(GriddingResult&)> finishCallback);

  /**
   * Block until all tasks have finished.
   */
  virtual void Finish(){};

  /**
   * Make the gridding task manager according to the settings.
   */
  static std::unique_ptr<GriddingTaskManager> Make(const Settings& settings);

 protected:
  Resources GetResources() const;

  /**
   * Run the provided task with the specified resources.
   * @param task A possibly compound gridding task.
   *        RunDirect() moves large values for the given facets out of the task.
   * @param facet_indices A list with the indices of the facets which should be
   *        gridded. The sequential GriddingTaskManager supplies all facets.
   *        The parallel ThreadedScheduler supplies a single facet and does
   *        multiple RunDirect calls for different facets in parallel.
   * @param resources Resources for creating the gridder.
   * @param result [out] Storage for gridding result.
   *        result.facets should have an entry for each facet.
   *        RunDirect() only updates the result for the given facet indices.
   *        When the first facet index is 0, it also updates the generic
   *        result part, which is equal for all facets.
   * @param result_mutex Protects concurrent accesses to result fields that
   *        are updated for each facet.
   */
  void RunDirect(GriddingTask& task, const std::vector<size_t>& facet_indices,
                 const Resources& resources, GriddingResult& result,
                 std::mutex& result_mutex);

 private:
  std::unique_ptr<MSGridderBase> ConstructGridder(
      const Resources& resources) const;

  /** Initializes 'gridder' with values that are equal for all facets. */
  void InitializeGridderForTask(MSGridderBase& gridder,
                                const GriddingTask& task);

  /** Initializes 'gridder' with facet-specific values. */
  void InitializeGridderForFacet(MSGridderBase& gridder,
                                 GriddingTask::FacetData& facet_task);

  const Settings& settings_;

  /**
   * Writer lock manager for the scheduler.
   * By default, it equals the 'this' pointer.
   * When the GriddingTaskManager is used within an MPIScheduler,
   * it may point to the MPIScheduler.
   */
  GriddingTaskManager* writer_lock_manager_;
};

#endif
