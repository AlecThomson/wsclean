#ifndef SCHEDULING_MPI_WORKER_SCHEDULER_H_
#define SCHEDULING_MPI_WORKER_SCHEDULER_H_

#include "griddingtaskmanager.h"

#include "../main/settings.h"

#include "griddingresult.h"

class MpiWorkerScheduler final : public GriddingTaskManager {
 public:
  MpiWorkerScheduler(const class Settings& settings);

  ~MpiWorkerScheduler() override { Finish(); }

  int Rank() const { return rank_; }

  /**
   * Run function for use in Worker.
   * Runs the task using the local scheduler.
   */
  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finish_callback) override {
    local_scheduler_->Run(std::move(task), finish_callback);
  }

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) override {
    return std::make_unique<WorkerWriterLock>(*this, writer_group_index);
  }

 private:
  class WorkerWriterLock : public WriterLock {
   public:
    explicit WorkerWriterLock(MpiWorkerScheduler& scheduler,
                              size_t writer_group_index);
    ~WorkerWriterLock() override;

   private:
    MpiWorkerScheduler&
        scheduler_;              ///< For accessing the scheduler's internals.
    size_t writer_group_index_;  ///< Index of the lock that must be acquired.
  };

  friend class WorkerWriterLock;

  /** MPI rank / node index. */
  int rank_;

  /**
   * The lower-level local scheduler on an MPI node.
   */
  std::unique_ptr<GriddingTaskManager> local_scheduler_;
};

#endif  // SCHEDULING_MPI_WORKER_SCHEDULER_H_
