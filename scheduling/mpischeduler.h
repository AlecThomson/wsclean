#ifndef SCHEDULING_MPI_SCHEDULER_H_
#define SCHEDULING_MPI_SCHEDULER_H_

#ifdef HAVE_MPI

#include "griddingtaskmanager.h"
#include "threadedscheduler.h"

#include <aocommon/queue.h>

#include <mutex>
#include <thread>
#include <condition_variable>

namespace wsclean {

class MPIScheduler final : public GriddingTaskManager {
 public:
  MPIScheduler(const class Settings& settings);
  ~MPIScheduler() { Finish(); }

  /**
   * Main Run function.
   * Either sends the task to another MPI node or runs it locally.
   */
  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finishCallback) override;

  void Finish() override;

  void Start(size_t nWriterGroups) override;

  std::unique_ptr<WriterLock> GetLock(size_t writerGroupIndex) override {
    // Since WSClean uses a static outputchannel-to-node mapping,
    // synchronisation of writes only needs to happen within a node.
    return _localScheduler.GetLock(writerGroupIndex);
  }

  /**
   * Send a task to a worker node or run it on the master
   * If all nodes are busy, the call will block until a node is available.
   */
  void send(GriddingTask&& task,
            std::function<void(GriddingResult&)>&& callback);

  /**
   * Wait until results are available and push these to the 'ready list'.
   * The loop ends when Finish() is called and all tasks are finished.
   * This function runs in a separate thread.
   */
  void receiveLoop();

  /**
   * Gets a node index for executing a (compound) task according
   * to the channel to index mapping.
   * Updates the number of tasks assigned to the node and stores
   * the callback function.
   * @return The index of the node executing the task.
   */
  int getNode(const GriddingTask& task,
              std::function<void(GriddingResult&)>&& callback);

  /**
   * If any results are available, call the callback functions and remove these
   * results from the ready list. This function
   * should be called by the main thread only, so that the user of the MPI
   * scheduler does not need to synchronize.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  void processReadyList_UNSYNCHRONIZED();

  /**
   * Return true if any tasks are still running.
   * Remember that the return value is independent of the state of the
   * master node: when the master node is gridding, it will nevertheless
   * return false if the other nodes are not running tasks.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  bool receiveTasksAreRunning_UNSYNCHRONIZED();

  void processGriddingResult(int node, size_t bodySize);
  /**
   * Stores 'result' in _readyList and updates the available slots of 'node'.
   */
  void StoreResult(GriddingResult&& result, int node);

  bool _isRunning;
  bool _isFinishing;
  std::condition_variable _notify;
  std::mutex _mutex;
  std::thread _receiveThread;
  /** Stores results of ready tasks. */
  std::vector<GriddingResult> _readyList;
  /** Stores callbacks, indexed by task id. */
  std::map<size_t, std::function<void(GriddingResult&)>> _callbacks;

  /**
   * Available execution room for tasks for each node.
   * A compound task, with multiple facets, counts as n_facets tasks.
   * The value is negative if the scheduler sent more tasks to the
   * node than it can execute in parallel. This situation occurs if:
   * - A compound task contains more sub-tasks than the available room.
   *   The scheduler still schedules such tasks, since there's no need to wait
   *   until the node has room for the entire task. In it's available room,
   *   the node can start with sub-tasks instead of being idle.
   * - The scheduler sends tasks prematurely. When a node finishes a task,
   *   it can then immediately start with the prematurely sent task instead
   *   of waiting for a new task.
   */
  std::vector<int> _availableRoom;

  /**
   * The lower-level local scheduler on an MPI node.
   * Using the threaded scheduler ensures that gridding uses a separate thread.
   */
  ThreadedScheduler _localScheduler;
};

#endif  // HAVE_MPI

}  // namespace wsclean

#endif  // MPI_SCHEDULER_H
