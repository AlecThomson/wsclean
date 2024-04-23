#include "mpischeduler.h"

#include <algorithm>
#include <cassert>
#include <memory>

#include "griddingresult.h"

#include "../main/settings.h"

#include "../distributed/mpibig.h"
#include "../distributed/taskmessage.h"

#include <aocommon/logger.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

using aocommon::Logger;

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
constexpr int kSlotsPerNode = 1;
}  // namespace

MPIScheduler::MPIScheduler(const Settings& settings)
    : GriddingTaskManager(settings),
      _isRunning(false),
      _isFinishing(false),
      _mutex(),
      _send_mutex(),
      _receiveThread(),
      _readyList(),
      _callbacks(),
      _availableRoom(),
      _writerLockQueues(),
      _localScheduler(settings) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  _availableRoom.assign(world_size, settings.parallelGridding);
  if (!settings.masterDoesWork) {
    _availableRoom[0] = 0;
  }
  _localScheduler.SetWriterLockManager(*this);
}

void MPIScheduler::Run(GriddingTask&& task,
                       std::function<void(GriddingResult&)> finishCallback) {
  if (!_isRunning) {
    _isFinishing = false;
    if (_availableRoom.size() > 1)
      _receiveThread = std::thread([&]() { receiveLoop(); });
    _isRunning = true;
  }
  send(std::move(task), std::move(finishCallback));

  std::lock_guard<std::mutex> lock(_mutex);
  processReadyList_UNSYNCHRONIZED();
}

void MPIScheduler::Finish() {
  if (_isRunning) {
    Logger::Info << "Finishing scheduler.\n";
    _localScheduler.Finish();

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

    if (_availableRoom.size() > 1) _receiveThread.join();

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
  GriddingTaskManager::Start(nWriterGroups);
  _writerLockQueues.resize(nWriterGroups);
}

std::unique_ptr<GriddingTaskManager::WriterLock> MPIScheduler::GetLock(
    size_t writer_group_index) {
  assert(writer_group_index < _writerLockQueues.size());
  return std::make_unique<MasterWriterLock>(*this, writer_group_index);
}

void MPIScheduler::send(GriddingTask&& task,
                        std::function<void(GriddingResult&)>&& callback) {
  int node = getNode(task, std::move(callback));
  Logger::Info << "Sending gridding task " << task.unique_id << " to node "
               << node << ".\n";

  if (node == 0) {
    _localScheduler.Run(std::move(task), [this](GriddingResult& result) {
      Logger::Info << "Main node has finished gridding task "
                   << result.unique_id << ".\n";
      StoreResult(std::move(result), 0);
    });
  } else {
    aocommon::SerialOStream payloadStream;
    // To use MPI_Send_Big, a uint64_t need to be reserved
    payloadStream.UInt64(0);
    task.Serialize(payloadStream);

    TaskMessage message;
    message.type = TaskMessage::Type::kGriddingRequest;
    message.bodySize = payloadStream.size();

    aocommon::SerialOStream taskMessageStream;
    message.Serialize(taskMessageStream);
    assert(taskMessageStream.size() == TaskMessage::kSerializedSize);

    std::unique_lock<std::mutex> lock(_send_mutex);
    MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
             0, MPI_COMM_WORLD);
    MPI_Send_Big(payloadStream.data(), payloadStream.size(), node, 0,
                 MPI_COMM_WORLD, GetSettings().maxMpiMessageSize);
  }
}

int MPIScheduler::getNode(const GriddingTask& task,
                          std::function<void(GriddingResult&)>&& callback) {
  std::unique_lock<std::mutex> lock(_mutex);

  assert(_callbacks.count(task.unique_id) == 0);
  _callbacks.emplace(task.unique_id, callback);

  while (true) {
    // Find the node that has the most available execution slots.
    // The backwards search prefers worker nodes over the main node.
    const auto nodeWithMostSlots =
        std::max_element(_availableRoom.rbegin(), _availableRoom.rend());
    // Select nodeWithMostSlots even if it can't process all facets in parallel.
    // It's still the best candidate. Also, processing all facets in parallel
    // may not be possible at all.
    if (*nodeWithMostSlots > 0) {
      *nodeWithMostSlots -= task.facets.size();
      _notify.notify_all();
      return _availableRoom.rend() - nodeWithMostSlots - 1;
    } else {
      // All nodes are busy -> Wait and try again later.
      _notify.wait(lock);
    }
  }
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

      const int node = status.MPI_SOURCE;
      switch (message.type) {
        case TaskMessage::Type::kGriddingResult:
          processGriddingResult(node, message.bodySize);
          break;
        case TaskMessage::Type::kLockRequest:
          processLockRequest(node, message.lockId);
          break;
        case TaskMessage::Type::kLockRelease:
          processLockRelease(node, message.lockId);
          break;
        default:
          throw std::runtime_error("Invalid message sent by node " +
                                   std::to_string(node));
      }

      lock.lock();
    }
  }
  Logger::Info << "All worker nodes have finished their gridding tasks.\n";
}

void MPIScheduler::processReadyList_UNSYNCHRONIZED() {
  while (!_readyList.empty()) {
    // Call the callback for this finished task
    GriddingResult& result = _readyList.back();
    // Copy the task id, since callbacks may adjust the result.
    const size_t task_id = result.unique_id;
    _callbacks[task_id](result);
    _readyList.pop_back();
    _callbacks.erase(task_id);
  }
}

bool MPIScheduler::receiveTasksAreRunning_UNSYNCHRONIZED() {
  for (size_t i = 1; i != _availableRoom.size(); ++i) {
    if (_availableRoom[i] < static_cast<int>(GetSettings().parallelGridding)) {
      return true;
    }
  }
  return false;
}

void MPIScheduler::processGriddingResult(int node, size_t bodySize) {
  aocommon::UVector<unsigned char> buffer(bodySize);
  MPI_Status status;
  MPI_Recv_Big(buffer.data(), bodySize, node, 0, MPI_COMM_WORLD, &status,
               GetSettings().maxMpiMessageSize);
  GriddingResult result;
  aocommon::SerialIStream stream(std::move(buffer));
  stream.UInt64();  // storage for MPI_Recv_Big
  result.Unserialize(stream);
  StoreResult(std::move(result), node);
}

void MPIScheduler::StoreResult(GriddingResult&& result, int node) {
  std::lock_guard<std::mutex> lock(_mutex);
  _availableRoom[node] += result.facets.size();
  _readyList.emplace_back(std::move(result));
  _notify.notify_all();
}

void MPIScheduler::processLockRequest(int node, size_t lockId) {
  if (lockId >= _writerLockQueues.size()) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " requests invalid lock " +
                             std::to_string(lockId));
  }

  std::unique_lock<std::mutex> lock(_mutex);

  // Idle nodes shouldn't request locks.
  assert(_availableRoom[node] <
         static_cast<int>(GetSettings().parallelGridding));

  _writerLockQueues[lockId].PushBack(node);

  if (_writerLockQueues[lockId].Size() == 1) {
    // Grant the lock immediately.
    lock.unlock();
    grantLock(node, lockId);
  }  // else processLockRelease will grant the lock once it's released.
}

void MPIScheduler::processLockRelease(int node, size_t lockId) {
  if (lockId >= _writerLockQueues.size()) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " releases invalid lock id " +
                             std::to_string(lockId));
  }

  std::unique_lock<std::mutex> lock(_mutex);
  if (_writerLockQueues[lockId].Empty() ||
      _writerLockQueues[lockId][0] != node) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " releases not-granted lock id " +
                             std::to_string(lockId));
  }

  _writerLockQueues[lockId].PopFront();
  if (!_writerLockQueues[lockId].Empty()) {
    const int waiting_node = _writerLockQueues[lockId][0];
    if (waiting_node == 0) {
      // Notify the worker thread, which waits in MasterWriterLock::lock.
      _notify.notify_all();
    } else {
      lock.unlock();
      grantLock(waiting_node, lockId);
    }
  }
}

void MPIScheduler::grantLock(int node, size_t lockId) {
  const TaskMessage message(TaskMessage::Type::kLockGrant, lockId);
  aocommon::SerialOStream taskMessageStream;
  message.Serialize(taskMessageStream);
  assert(taskMessageStream.size() == TaskMessage::kSerializedSize);

  // Using asynchronous MPI_ISend is possible here, however, a synchronous
  // MPI_Send is much simpler. Since the message is small and the receiver
  // is already waiting for the message, the overhead of synchronous
  // communication should be limited.
  std::unique_lock<std::mutex> lock(_send_mutex);
  MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
           0, MPI_COMM_WORLD);
}

MPIScheduler::MasterWriterLock::MasterWriterLock(MPIScheduler& scheduler,
                                                 size_t writer_group_index)
    : scheduler_(scheduler), writer_group_index_(writer_group_index) {
  std::unique_lock<std::mutex> lock(scheduler_._mutex);
  assert(scheduler_._availableRoom[0] <
         static_cast<int>(scheduler_.GetSettings().parallelGridding));
  scheduler_._writerLockQueues[writer_group_index].PushBack(0);

  while (scheduler_._writerLockQueues[writer_group_index][0] != 0) {
    // Wait for lock in the worker thread on the master node.
    // processLockRelease notifies the worker thread when the lock is available.
    scheduler_._notify.wait(lock);
  }
}

MPIScheduler::MasterWriterLock::~MasterWriterLock() {
  scheduler_.processLockRelease(0, writer_group_index_);
}
