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
  const std::size_t task_id = task.unique_id;
  const std::size_t facet_count = task.facets.size();
  TaskData* task_data;

  // Add an entry in task_data_map_ for the task.
  {
    std::lock_guard<std::mutex> lock(mutex_);
    assert(task_data_map_.count(task_id) == 0);
    task_data = &task_data_map_[task_id];
  }
  task_data->task = std::move(task);
  task_data->result.facets.resize(facet_count);
  task_data->callback = std::move(finish_callback);

  // Add sub-tasks for each facet to the task queue.
  for (std::size_t facet_index = 0; facet_index < facet_count; ++facet_index) {
    task_queue_.Emplace(task_id, facet_index);
  }

  ProcessReadyList();
}

void ThreadedScheduler::ProcessQueue() {
  std::pair<size_t, size_t> facet_task_pair;
  while (task_queue_.Pop(facet_task_pair)) {
    const std::size_t task_id = facet_task_pair.first;
    const std::size_t facet_index = facet_task_pair.second;
    TaskData& task_data = task_data_map_[task_id];
    try {
      RunDirect(task_data.task, {facet_index}, resources_per_task_,
                task_data.result, task_data.result_mutex);
    } catch (std::exception&) {
      std::lock_guard<std::mutex> lock(mutex_);
      latest_exception_ = std::current_exception();
    }

    // Extract the new value from task_data.finished_facet_count directly.
    // When extracting the new value later, another thread may also have
    // incremented the atomic value in the mean time, and multiple threads
    // will think they have processed the last facet.
    const size_t finished_facet_count = ++task_data.finished_facet_count;
    if (finished_facet_count == task_data.task.facets.size()) {
      if (GetSettings().UseMpi()) {
        // Execute callback immediately, from the processing thread.
        // The MPIScheduler stores the result at the main node.
        // The MPIWorkerScheduler sends the result to the main node.
        task_data.callback(task_data.result);
        std::lock_guard<std::mutex> lock(mutex_);
        task_data_map_.erase(task_id);
      } else {
        // Store the task id and execute the callback on the main thread.
        std::lock_guard<std::mutex> lock(mutex_);
        ready_list_.emplace_back(task_id);
      }
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
  std::vector<std::size_t> local_ready_list;
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
  for (std::size_t task_id : local_ready_list) {
    TaskData& task_data = task_data_map_[task_id];
    task_data.callback(task_data.result);
  }

  // Remove the finished tasks from task_data_map_.
  if (!local_ready_list.empty()) {
    std::lock_guard<std::mutex> lock{mutex_};
    for (std::size_t task_id : local_ready_list) {
      task_data_map_.erase(task_id);
    }
  }

  if (local_exception) std::rethrow_exception(local_exception);
}
