#include "mpiworkerscheduler.h"

#include "../distributed/taskmessage.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
}  // namespace

MpiWorkerScheduler::MpiWorkerScheduler(const Settings& settings)
    : GriddingTaskManager(settings),
      rank_(-1),
      local_scheduler_(GriddingTaskManager::MakeLocal(settings)) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  local_scheduler_->SetWriterLockManager(*this);
}

MpiWorkerScheduler::WorkerWriterLock::WorkerWriterLock(
    MpiWorkerScheduler& scheduler, size_t writer_group_index)
    : scheduler_(scheduler), writer_group_index_(writer_group_index) {
  TaskMessage message(TaskMessage::Type::kLockRequest, writer_group_index);
  aocommon::SerialOStream task_message_stream;
  message.Serialize(task_message_stream);

  MPI_Send(task_message_stream.data(), task_message_stream.size(), MPI_BYTE,
           kMainNode, kTag, MPI_COMM_WORLD);

  MPI_Status status;
  aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
  MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE, kMainNode,
           kTag, MPI_COMM_WORLD, &status);
  aocommon::SerialIStream stream(std::move(buffer));
  message.Unserialize(stream);
  if (message.type != TaskMessage::Type::kLockGrant ||
      message.lockId != writer_group_index) {
    int node = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    throw std::runtime_error("Node " + std::to_string(node) +
                             " received an invalid message from node " +
                             std::to_string(status.MPI_SOURCE) + " (type " +
                             std::to_string(static_cast<int>(message.type)) +
                             ", lock " + std::to_string(message.lockId) +
                             ") while requesting lock " +
                             std::to_string(writer_group_index));
  }
}

MpiWorkerScheduler::WorkerWriterLock::~WorkerWriterLock() {
  TaskMessage message(TaskMessage::Type::kLockRelease, writer_group_index_);
  aocommon::SerialOStream task_message_stream;
  message.Serialize(task_message_stream);

  // Using asynchronous MPI_ISend is possible here, however, the
  // task_message_stream should then remain valid after the call and the extra
  // 'request' handle probably needs handling.
  MPI_Send(task_message_stream.data(), task_message_stream.size(), MPI_BYTE,
           kMainNode, kTag, MPI_COMM_WORLD);
}
