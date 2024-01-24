#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include "griddingtask.h"
#include "griddingresult.h"

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"
#include "../structures/msselection.h"

#include <aocommon/polarization.h>

#include "../msproviders/msdatadescription.h"

#include <aocommon/lane.h>

#include <cstring>
#include <functional>
#include <memory>
#include <vector>

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

  void SetWriterLockManager(GriddingTaskManager& manager) {
    _writerLockManager = &manager;
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

  const Settings& _settings;

  std::unique_ptr<MSGridderBase> makeGridder(const Resources& resources) const;

  /**
   * Run the provided task with the specified gridder.
   */
  GriddingResult RunDirect(GriddingTask&& task, MSGridderBase& gridder);

  /**
   * Make a local gridding task manager according to the settings.
   * Both Make() and the MPIScheduler use this function. The MPIScheduler uses
   * a local scheduler for running tasks at an MPI node.
   */
  static std::unique_ptr<GriddingTaskManager> MakeLocal(
      const Settings& settings);

 private:
  std::unique_ptr<MSGridderBase> constructGridder(
      const Resources& resources) const;

  /**
   * Writer lock manager for the scheduler.
   * By default, it equals the 'this' pointer.
   * When the GriddingTaskManager is used within an MPIScheduler,
   * it may point to the MPIScheduler.
   */
  GriddingTaskManager* _writerLockManager;
};

#endif
