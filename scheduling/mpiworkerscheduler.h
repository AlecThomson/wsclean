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

  /**
   * Run function for use in Worker.
   * Runs the task using the local scheduler.
   * Note: MpiWorkerScheduler ignores the callback function.
   */
  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> ignored_callback) override;

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) override {
    return std::make_unique<WorkerWriterLock>(*this, writer_group_index);
  }

  void GrantLock(size_t writer_group_index);

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

  /** Common part of GrantLock() and ~WorkerWriterLock().*/
  void GrantLockUnsynchronized(size_t writer_group_index);

  /** MPI rank / node index. */
  int rank_;

  /**
   * Serializes MPI_Send calls from different threads.
   * Protects all *writer_locks_ members.
   */
  std::mutex mutex_;
  /** Synchronizes threads waiting for writer locks. */
  std::condition_variable notify_;
  /**
   * Signals that a writer lock is available for a thread.
   * A writer lock is in the set after receiving a lock-grant message or when
   * a thread released the lock while another thread is waiting for it.
   */
  std::set<size_t> available_writer_locks_;
  /**
   * For each writer lock, the number of local threads that either have or are
   * waiting for the writer lock.
   */
  std::map<size_t, size_t> writer_locks_;

  /**
   * The lower-level local scheduler on an MPI node.
   * Always use a ThreadedScheduler since acquiring writer locks in the gridder
   * should use a different thread than Worker::Run().
   */
  ThreadedScheduler local_scheduler_;
};

#endif  // SCHEDULING_MPI_WORKER_SCHEDULER_H_
