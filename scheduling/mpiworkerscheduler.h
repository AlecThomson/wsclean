#ifndef SCHEDULING_MPI_WORKER_SCHEDULER_H_
#define SCHEDULING_MPI_WORKER_SCHEDULER_H_

#include <condition_variable>
#include <mutex>
#include <set>

#include "griddingtaskmanager.h"
#include "threadedscheduler.h"

#include "../main/settings.h"

#include "griddingresult.h"

class MpiWorkerScheduler final : public GriddingTaskManager {
 public:
  MpiWorkerScheduler(const class Settings& settings);

  ~MpiWorkerScheduler() override { Finish(); }

  int Rank() const { return rank_; }

  void Start(size_t n_writer_groups) override {
    GriddingTaskManager::Start(n_writer_groups);
    local_scheduler_.Start(n_writer_groups);
  }

  /**
   * Run function for use in Worker.
   * Runs the task using the local scheduler.
   * Note: MpiWorkerScheduler ignores the callback function.
   */
  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> ignored_callback) override;

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) override {
    // Since WSClean uses a static outputchannel-to-node mapping,
    // synchronisation of writes only needs to happen within a node.
    return local_scheduler_.GetLock(writer_group_index);
  }

  void GrantLock(size_t writer_group_index);

 private:
  /** MPI rank / node index. */
  int rank_;

  /** Serializes MPI_Send calls from different threads. */
  std::mutex mutex_;

  /**
   * The lower-level local scheduler on an MPI node.
   * Always use a ThreadedScheduler since acquiring writer locks in the gridder
   * should use a different thread than Worker::Run().
   */
  ThreadedScheduler local_scheduler_;
};

#endif  // SCHEDULING_MPI_WORKER_SCHEDULER_H_
