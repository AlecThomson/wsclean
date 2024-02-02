#include "threadedscheduler.h"
#include "../gridding/msgridderbase.h"

#include "../main/settings.h"

#include <aocommon/logger.h>

#include <string>

ThreadedScheduler::ThreadedScheduler(const Settings& settings)
    : GriddingTaskManager{settings},
      // When using the ThreadedScheduler as the main scheduler, limit the
      // number of tasks in the queue to one per thread. When stacking too many
      // tasks, memory usage could become an issue.
      // When using the ThreadedScheduler as a local scheduler with the
      // MPIScheduler, the MPIScheduler manages the task distribution.
      // The ThreadedScheduler should always queue new tasks in that case.
      task_queue_(settings.UseMpi() ? TaskQueueType()
                                    : TaskQueueType(settings.parallelGridding)),
      resources_per_task_(GetResources().GetPart(settings.parallelGridding)) {
  for (size_t i = 0; i < settings.parallelGridding; ++i) {
    thread_list_.emplace_back(&ThreadedScheduler::ProcessQueue, this);
  }
}

ThreadedScheduler::~ThreadedScheduler() {
  try {
    Finish();
  } catch (std::exception& e) {
    // Normally, the user of the ThreadedScheduler calls Finish(), which
    // rethrows any exception caught in a thread.
    // We are in a destructor, so all that can be done is report the error.
    using namespace std::string_literals;
    aocommon::Logger::Error
        << "Exception caught during destruction of ThreadedScheduler:\n"s +
               e.what() + '\n';
  }
  task_queue_.Finish();  // Make all threads exit.
  for (std::thread& thread : thread_list_) thread.join();
}

void ThreadedScheduler::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finish_callback) {
  task_queue_.Emplace(std::move(task), std::move(finish_callback));

  ProcessReadyList();
}

void ThreadedScheduler::ProcessQueue() {
  std::pair<GriddingTask, std::function<void(GriddingResult&)>> task_pair;
  while (task_queue_.Pop(task_pair)) {
    GriddingResult result;
    try {
      std::unique_ptr<MSGridderBase> gridder(makeGridder(resources_per_task_));
      result = RunDirect(std::move(task_pair.first), *gridder);
    } catch (std::exception&) {
      std::lock_guard<std::mutex> lock(mutex_);
      latest_exception_ = std::current_exception();
    }

    std::function<void(GriddingResult&)> callback = std::move(task_pair.second);
    if (GetSettings().UseMpi()) {
      // Execute callback immediately, from the processing thread.
      // The MPIScheduler stores the result at the main node.
      // The MPIWorkerScheduler sends the result to the main node.
      callback(result);
    } else {
      // Store the result and execute the callback on the main thread.
      std::lock_guard<std::mutex> lock(mutex_);
      ready_list_.emplace_back(std::move(result), std::move(callback));
    }
  }
}

void ThreadedScheduler::Start(size_t nWriterGroups) {
  assert(ready_list_.empty());
  GriddingTaskManager::Start(nWriterGroups);
  if (writer_group_locks_.size() < nWriterGroups)
    writer_group_locks_ = std::vector<std::mutex>(nWriterGroups);
}

std::unique_ptr<GriddingTaskManager::WriterLock> ThreadedScheduler::GetLock(
    size_t writer_group_index) {
  assert(writer_group_index < writer_group_locks_.size());
  return std::make_unique<ThreadedWriterLock>(*this, writer_group_index);
}

void ThreadedScheduler::Finish() {
  task_queue_.WaitForIdle(GetSettings().parallelGridding);

  ProcessReadyList();
}

void ThreadedScheduler::ProcessReadyList() {
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult&)>>>
      local_ready_list;
  std::exception_ptr local_exception;

  // Move the ready_list_ and latest_exception_ to local variables so we can
  // release the lock while the callbacks run / while throwing the exception.
  {
    std::lock_guard<std::mutex> lock{mutex_};

    local_ready_list = std::move(ready_list_);
    ready_list_.clear();

    local_exception = std::move(latest_exception_);
    latest_exception_ = std::exception_ptr();
  }

  // Call callbacks for any finished tasks
  for (std::pair<GriddingResult, std::function<void(GriddingResult&)>>& item :
       local_ready_list) {
    item.second(item.first);
  }

  if (local_exception) std::rethrow_exception(local_exception);
}
