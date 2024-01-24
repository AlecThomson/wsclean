#ifndef SCHEDULING_THREADED_SCHEDULER_H_
#define SCHEDULING_THREADED_SCHEDULER_H_

#include <mutex>
#include <thread>

#include "griddingtaskmanager.h"

#include "../structures/resources.h"

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
  void CheckExceptions();

  std::mutex mutex_;
  std::exception_ptr latest_exception_;
  std::vector<std::thread> thread_list_;
  aocommon::Lane<std::pair<GriddingTask, std::function<void(GriddingResult&)>>>
      task_list_;
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult&)>>>
      ready_list_;
  std::vector<std::mutex> writer_group_locks_;

  const Resources resources_per_task_;
};

#endif
