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
      _receiveThread(),
      _readyList(),
      _callbacks(),
      _availableRoom(),
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

  const uint32_t mpi_writer_groups = nWriterGroups;
  for (size_t rank = 1; rank < _availableRoom.size(); ++rank) {
    assert(_availableRoom[rank] > 0 &&
           size_t(_availableRoom[rank]) == GetSettings().parallelGridding);
    MPI_Send(&mpi_writer_groups, 1, MPI_UINT32_T, rank, kTag, MPI_COMM_WORLD);
  }

  if (GetSettings().masterDoesWork) {
    _localScheduler.Start(nWriterGroups);
  }
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

    MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
             0, MPI_COMM_WORLD);
    MPI_Send_Big(payloadStream.data(), payloadStream.size(), node, 0,
                 MPI_COMM_WORLD, GetSettings().maxMpiMessageSize);
  }
}

int MPIScheduler::getNode(const GriddingTask& task,
                          std::function<void(GriddingResult&)>&& callback) {
  // Determine the target node using the channel to node mapping.
  int node = GetSettings().channelToNode[task.outputChannelIndex];

  // Wait until _availableRoom[node] becomes larger than 0.
  std::unique_lock<std::mutex> lock(_mutex);
  while (_availableRoom[node] <= 0) {
    _notify.wait(lock);
  }
  _availableRoom[node] -= task.facets.size();

  // Store the callback function.
  assert(_callbacks.count(task.unique_id) == 0);
  _callbacks.emplace(task.unique_id, std::move(callback));

  return node;
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
