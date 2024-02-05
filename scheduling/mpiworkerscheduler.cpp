#include "mpiworkerscheduler.h"

#include <mpi.h>

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>
#include <aocommon/logger.h>

#include "../distributed/taskmessage.h"
#include "../distributed/mpibig.h"

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
}  // namespace

MpiWorkerScheduler::MpiWorkerScheduler(const Settings& settings)
    : GriddingTaskManager{settings}, rank_{-1}, local_scheduler_{settings} {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  local_scheduler_.SetWriterLockManager(*this);
}

void MpiWorkerScheduler::Run(
    GriddingTask&& task,
    [[maybe_unused]] std::function<void(GriddingResult&)> ignored_callback) {
  aocommon::Logger::Info << "Worker node " << rank_
                         << " is starting gridding.\n";
  local_scheduler_.Run(std::move(task), [this](GriddingResult& result) {
    aocommon::Logger::Info << "Worker node " << rank_ << " is done gridding.\n";

    aocommon::SerialOStream resStream;
    resStream.UInt64(0);  // reserve nr of packages for MPI_Send_Big
    result.Serialize(resStream);

    TaskMessage message;
    message.type = TaskMessage::Type::kGriddingResult;
    message.bodySize = resStream.size();

    aocommon::SerialOStream msgStream;
    message.Serialize(msgStream);
    assert(msgStream.size() == TaskMessage::kSerializedSize);

    std::lock_guard<std::mutex> lock(mutex_);
    MPI_Send(msgStream.data(), msgStream.size(), MPI_BYTE, 0, 0,
             MPI_COMM_WORLD);
    MPI_Send_Big(resStream.data(), resStream.size(), 0, 0, MPI_COMM_WORLD);
  });
}

void MpiWorkerScheduler::GrantLock(size_t writer_group_index) {
  std::lock_guard<std::mutex> lock{mutex_};
  assert(writer_locks_.count(writer_group_index) == 0);
  writer_locks_.insert(writer_group_index);
  notify_.notify_all();
}

MpiWorkerScheduler::WorkerWriterLock::WorkerWriterLock(
    MpiWorkerScheduler& scheduler, size_t writer_group_index)
    : scheduler_(scheduler), writer_group_index_(writer_group_index) {
  TaskMessage message(TaskMessage::Type::kLockRequest, writer_group_index);
  aocommon::SerialOStream task_message_stream;
  message.Serialize(task_message_stream);

  std::unique_lock<std::mutex> lock{scheduler.mutex_};
  assert(scheduler_.writer_locks_.count(writer_group_index) == 0);
  MPI_Send(task_message_stream.data(), task_message_stream.size(), MPI_BYTE,
           kMainNode, kTag, MPI_COMM_WORLD);
  do {
    scheduler_.notify_.wait(lock);
  } while (scheduler_.writer_locks_.count(writer_group_index) == 0);
}

MpiWorkerScheduler::WorkerWriterLock::~WorkerWriterLock() {
  TaskMessage message(TaskMessage::Type::kLockRelease, writer_group_index_);
  aocommon::SerialOStream task_message_stream;
  message.Serialize(task_message_stream);

  std::lock_guard<std::mutex> lock{scheduler_.mutex_};
  // Using asynchronous MPI_ISend is possible here, however, the
  // task_message_stream should then remain valid after the call and the extra
  // 'request' handle probably needs handling.
  MPI_Send(task_message_stream.data(), task_message_stream.size(), MPI_BYTE,
           kMainNode, kTag, MPI_COMM_WORLD);
  [[maybe_unused]] const size_t erased =
      scheduler_.writer_locks_.erase(writer_group_index_);
  assert(erased == 1);
}
