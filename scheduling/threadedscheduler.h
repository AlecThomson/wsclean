#ifndef SCHEDULING_THREADED_SCHEDULER_H_
#define SCHEDULING_THREADED_SCHEDULER_H_

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>

#include <aocommon/taskqueue.h>

#include "griddingtaskmanager.h"

#include "../structures/resources.h"

namespace wsclean {

class ThreadedScheduler final : public GriddingTaskManager {
 public:
  ThreadedScheduler(const class Settings& settings);
  ~ThreadedScheduler();

  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finishCallback) override;
  void Finish() override;

  void Start(size_t nWriterGroups) override;

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) override;

 private:
  class ThreadedWriterLock final : public WriterLock {
   public:
    explicit ThreadedWriterLock(ThreadedScheduler& scheduler,
                                size_t writer_group_index)
        : scheduler_{scheduler}, writer_group_index_{writer_group_index} {
      scheduler.writer_group_locks_[writer_group_index].lock();
    }

    ~ThreadedWriterLock() override {
      scheduler_.writer_group_locks_[writer_group_index_].unlock();
    }

   private:
    ThreadedScheduler& scheduler_;
    size_t writer_group_index_;
  };

  friend class ThreadedWriterLock;

  void ProcessQueue();
  void ProcessReadyList();

  /// Contains all data for a single task.
  struct TaskData {
    GriddingTask task;
    GriddingResult result;
    std::mutex result_mutex;
    std::atomic<std::size_t> finished_facet_count;
    std::function<void(GriddingResult&)> callback;
  };

  /// The first value of the pair is the unique id for the tasks, for looking up
  /// the task and its callback into task_data_map_.
  /// The second value is the index of the facet for the sub-task.
  using TaskQueueType = aocommon::TaskQueue<std::pair<size_t, size_t>>;

  /// Protects task_data_map_, latest_exception_ and ready_list_.
  std::mutex mutex_;
  std::vector<std::thread> thread_list_;
  /// Stores the data for the pending tasks.
  /// The map is indexed by the unique_id of the task.
  std::map<std::size_t, TaskData> task_data_map_;
  /// FIFO queue for the tasks.
  /// Each sub-task, which processes a single facet, has an entry.
  TaskQueueType task_queue_;
  /// Stores the unique_id for tasks that are complete, so the main thread
  /// can execute the callback.
  std::vector<size_t> ready_list_;
  /// Stores the latest exception that occurred while running a task in a
  /// thread, so the scheduler can rethrow it in the main thread.
  std::exception_ptr latest_exception_;
  std::vector<std::mutex> writer_group_locks_;

  const Resources resources_per_task_;
};

}  // namespace wsclean

#endif
