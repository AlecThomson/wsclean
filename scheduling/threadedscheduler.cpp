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

  // Set the amount of parallel gridders for this task
  task.num_parallel_gridders_ = std::min(
      GetSettings().parallelGridding, std::max(task.facets.size(), size_t{1}));

  // Add an entry in task_data_map_ for the task.
  {
    std::lock_guard<std::mutex> lock(mutex_);
    assert(task_data_map_.count(task_id) == 0);
    task_data = &task_data_map_[task_id];
  }
  task_data->task = std::move(task);
  task_data->result.facets.resize(facet_count);
  task_data->callback = std::move(finish_callback);

  if (task_data->task.operation == GriddingTask::Invert) {
    std::lock_guard<std::mutex> lock(mutex_);

    // There must be N-1 blocker tasks before we can run, giving us N threads in
    // total as requested
    // Set up locks so that we only start operating once we have the resources
    // The scheduler will allocate us resources when the dummy tasks are placed
    // in the queue
    task_data->task.batch_task_mutex = std::make_shared<std::mutex>();
    task_data->task.batch_task_mutex->lock();

    for (size_t block_task_index = 0;
         block_task_index < task_data->task.num_parallel_gridders_ - 1;
         ++block_task_index) {
      size_t block_task_id = std::numeric_limits<size_t>::max() -
                             (task_id * 1000) - block_task_index;
      TaskData* block_task_data;

      assert(task_data_map_.count(block_task_id) == 0);
      block_task_data = &task_data_map_[block_task_id];

      GriddingTask block_task;
      block_task.operation = GriddingTask::Block;
      block_task.batch_task_mutex = task_data->task.batch_task_mutex;
      block_task_data->task = std::move(block_task);
    }
  }

  task_queue_.Emplace(task_id, std::numeric_limits<size_t>::max());

  if (task_data->task.operation == GriddingTask::Invert) {
    // Create some blocking tasks that will take a resource from the scheduler
    // (when available) and assign it back to us
    for (size_t block_task_index = 0;
         block_task_index < task_data->task.num_parallel_gridders_ - 1;
         ++block_task_index) {
      size_t block_task_id = std::numeric_limits<size_t>::max() -
                             (task_data->task.unique_id * 1000) -
                             block_task_index;
      task_queue_.Emplace(block_task_id, 0);
    }
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
      if (facet_index == std::numeric_limits<size_t>::max()) {
        std::vector<size_t> facet_indexes;
        for (const auto& facet : task_data.task.facets) {
          facet_indexes.push_back(facet.index);
        }
        Resources task_resources = resources_per_task_;
        task_resources =
            resources_per_task_.GetCombined(GetSettings().parallelGridding);
        RunDirect(task_data.task, facet_indexes, task_resources,
                  task_data.result, task_data.result_mutex);
      } else {
        RunDirect(task_data.task, {facet_index}, resources_per_task_,
                  task_data.result, task_data.result_mutex);
      }
    } catch (std::exception&) {
      std::lock_guard<std::mutex> lock(mutex_);
      latest_exception_ = std::current_exception();
    }

    if (task_data.task.operation == GriddingTask::Block) {
      std::lock_guard<std::mutex> lock(mutex_);
      task_data_map_.erase(task_id);
      continue;
    }

    if (facet_index == std::numeric_limits<size_t>::max()) {
      task_data.finished_facet_count += task_data.task.facets.size();
    } else {
      task_data.finished_facet_count++;
    }

    // Extract the new value from task_data.finished_facet_count directly.
    // When extracting the new value later, another thread may also have
    // incremented the atomic value in the mean time, and multiple threads
    // will think they have processed the last facet.
    bool process = false;
    task_data.result_mutex.lock();
    if (task_data.finished_facet_count == task_data.task.facets.size()) {
      process = true;
      task_data.finished_facet_count = 0;
    }
    task_data.result_mutex.unlock();
    if (process) {
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
