#include "mpischeduler.h"

#include "griddingresult.h"

#include "../io/logger.h"

#include "../main/settings.h"

#include "../distributed/mpibig.h"
#include "../distributed/taskmessage.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

#include <cassert>

MPIScheduler::MPIScheduler(const Settings &settings)
    : GriddingTaskManager(settings),
      _masterDoesWork(settings.masterDoesWork),
      _isRunning(false) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  _nodes.assign(
      world_size,
      std::make_pair(AvailableNode, std::function<void(GriddingResult &)>()));
  if (!settings.masterDoesWork && world_size <= 1)
    throw std::runtime_error(
        "Master was told not to work, but no other workers available");
}

MPIScheduler::~MPIScheduler() { Finish(); }

void MPIScheduler::Run(GriddingTask &&task,
                       std::function<void(GriddingResult &)> finishCallback) {
  if (!_isRunning) {
    _isFinishing = false;
    if (_nodes.size() > 1)
      _receiveThread = std::thread([&]() { receiveLoop(); });
    _isRunning = true;
  }
  send(std::move(task), finishCallback);

  std::lock_guard<std::mutex> lock(_mutex);
  processReadyList_UNSYNCHRONIZED();
}

void MPIScheduler::Finish() {
  Logger::Info << "Finishing scheduler.\n";
  if (_isRunning) {
    std::unique_lock<std::mutex> lock(_mutex);
    _isFinishing = true;
    _notify.notify_all();

    // As long as receive tasks are running, wait and keep processing
    // the ready list
    processReadyList_UNSYNCHRONIZED();
    while (receiveTasksAreRunning_UNSYNCHRONIZED()) {
      _notify.wait(lock);
      processReadyList_UNSYNCHRONIZED();
    }

    lock.unlock();

    if (_nodes.size() > 1) _receiveThread.join();

    if (_workThread.joinable()) _workThread.join();
    _isRunning = false;

    // The while loop above ignores the work thread, which might
    // be gridding on the master node. Therefore, the master thread
    // might have added an item to the ready list. Therefore,
    // the ready list should once more be processed.
    // A lock is no longer required, because all threads have stopped.
    processReadyList_UNSYNCHRONIZED();
  }
}

void MPIScheduler::Start(size_t nWriterGroups) {
  if (_nodes.size() > 1) {
    // FIXME: exception breaks also "non-facetting" mpiruns
    std::cout
        << "Requested an MPI lock for a run with #nodes > 1. This is not yet "
           "implemented"
        << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  GriddingTaskManager::Start(nWriterGroups);
  if (_writerGroupLocks.size() < nWriterGroups)
    _writerGroupLocks = std::vector<MPIWriterLock>(nWriterGroups);
}

WriterLockManager::LockGuard MPIScheduler::GetLock(size_t writerGroupIndex) {
  return LockGuard(_writerGroupLocks[writerGroupIndex]);
}

void MPIScheduler::runTaskOnNode0(GriddingTask &&task) {
  GriddingResult result = RunDirect(std::move(task));
  Logger::Info << "Master node has finished a gridding task.\n";
  std::unique_lock<std::mutex> lock(_mutex);
  _readyList.emplace_back(std::move(result), _nodes[0].second);
  _nodes[0].first = AvailableNode;
  lock.unlock();
  _notify.notify_all();
}

void MPIScheduler::send(GriddingTask &&task,
                        const std::function<void(GriddingResult &)> &callback) {
  int node =
      findAndSetNodeState(AvailableNode, std::make_pair(BusyNode, callback));
  Logger::Info << "Sending gridding task to node : " << node << '\n';

  if (node == 0) {
    if (_workThread.joinable()) _workThread.join();
    _workThread =
        std::thread(&MPIScheduler::runTaskOnNode0, this, std::move(task));
  } else {
    aocommon::SerialOStream payloadStream;
    // To use MPI_Send_Big, a uint64_t need to be reserved
    payloadStream.UInt64(0);
    task.Serialize(payloadStream);

    TaskMessage message;
    message.type = TaskMessage::GriddingRequest;
    message.bodySize = payloadStream.size();

    aocommon::SerialOStream taskMessageStream;
    message.Serialize(taskMessageStream);
    assert(taskMessageStream.size() == TaskMessage::kSerializedSize);

    MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
             0, MPI_COMM_WORLD);
    MPI_Send_Big(payloadStream.data(), payloadStream.size(), node, 0,
                 MPI_COMM_WORLD);
  }
}

int MPIScheduler::findAndSetNodeState(
    MPIScheduler::NodeState currentState,
    std::pair<MPIScheduler::NodeState, std::function<void(GriddingResult &)>>
        newState) {
  std::unique_lock<std::mutex> lock(_mutex);
  do {
    size_t iterEnd = _masterDoesWork ? _nodes.size() : _nodes.size() - 1;
    for (size_t i = 0; i != iterEnd; ++i) {
      const int node = _nodes.size() - i - 1;
      if (_nodes[node].first == currentState) {
        _nodes[node] = newState;
        _notify.notify_all();
        return node;
      }
    }
    _notify.wait(lock);
  } while (true);
}

void MPIScheduler::receiveLoop() {
  std::unique_lock<std::mutex> lock(_mutex);
  while (!_isFinishing || receiveTasksAreRunning_UNSYNCHRONIZED()) {
    if (!receiveTasksAreRunning_UNSYNCHRONIZED()) {
      _notify.wait(lock);
    } else {
      lock.unlock();

      TaskMessage message;
      MPI_Status status;
      aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
      MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE,
               MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      aocommon::SerialIStream stream(std::move(buffer));
      message.Unserialize(stream);

      int node = status.MPI_SOURCE;
      if (message.type != TaskMessage::GriddingResult)
        throw std::runtime_error("Invalid message sent by node " +
                                 std::to_string(node));

      buffer.resize(message.bodySize);
      MPI_Recv_Big(buffer.data(), message.bodySize, node, 0, MPI_COMM_WORLD,
                   &status);
      GriddingResult result;
      stream = aocommon::SerialIStream(std::move(buffer));
      stream.UInt64();  // storage for MPI_Recv_Big
      result.Unserialize(stream);

      lock.lock();
      _readyList.emplace_back(std::move(result), _nodes[node].second);
      _nodes[node].first = AvailableNode;
      lock.unlock();

      _notify.notify_all();

      lock.lock();
    }
  }
  Logger::Info << "All worker nodes have finished their gridding tasks.\n";
}

void MPIScheduler::processReadyList_UNSYNCHRONIZED() {
  while (!_readyList.empty()) {
    // Call the callback for this finished task
    _readyList.back().second(_readyList.back().first);
    _readyList.pop_back();
  }
}

bool MPIScheduler::receiveTasksAreRunning_UNSYNCHRONIZED() {
  for (size_t i = 1; i != _nodes.size(); ++i)
    if (_nodes[i].first == BusyNode) return true;
  return false;
}
