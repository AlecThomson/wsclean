#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include "griddingtask.h"
#include "griddingresult.h"
#include "griddinglockmanager.h"

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"
#include "../structures/msselection.h"

#include <aocommon/polarization.h>

#include "../msproviders/msdatadescription.h"

#include <aocommon/lane.h>

#include <cstring>
#include <functional>
#include <vector>

class MSGridderBase;

class GriddingTaskManager : public GriddingLockManager {
 public:
  virtual ~GriddingTaskManager();

  /**
   * Initialize writer groups. Call this function before scheduling Predict
   * tasks.
   *
   * This function resets the counter for all write groups to zero.
   * @param nWriterGroups The number of writer groups.
   */
  virtual void Start(size_t nWriterGroups) {
    _writerGroupCounters.assign(nWriterGroups, 0);
  }

  /**
   * Get the counter for a writer group and increment it for the next call.
   * Each call thus yields a higher value.
   * Note: Parallel GriddingTaskManager implementations should use proper
   * locking before using this base implementation.
   */
  WriterGroupLockGuard LockWriterGroup(size_t writerGroupIndex) const override {
    static DummyWriterLock dummy;
    WriterGroupLockGuard guard(dummy, _writerGroupCounters[writerGroupIndex]);
    return guard;
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
   * Run the given task. This variant of Run() does not call a
   * callback function when finished. A gridder is created for
   * the duration of the call.
   */
  GriddingResult RunDirect(GriddingTask&& task);

  /**
   * Make the gridding task manager according to the settings.
   */
  static std::unique_ptr<GriddingTaskManager> Make(
      const class Settings& settings, bool useDirectScheduler = false);

 protected:
  GriddingTaskManager(const class Settings& settings);

  const class Settings& _settings;

  size_t& getWriterGroupCounter(size_t writerGroupIndex) const {
    return _writerGroupCounters[writerGroupIndex];
  }

  std::unique_ptr<MSGridderBase> makeGridder() const;

  /**
   * Run the provided task with the specified gridder.
   */
  GriddingResult runDirect(GriddingTask&& task, MSGridderBase& gridder);

 private:
  class DummyWriterLock : public WriterLockBase {
   public:
    void lock() override final {}
    void unlock() override final {}
  };

  std::unique_ptr<MSGridderBase> constructGridder() const;

  mutable std::vector<size_t> _writerGroupCounters;
};

#endif
