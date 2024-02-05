#include "worker.h"

#include <mpi.h>

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/logger.h>

#include "../scheduling/griddingtask.h"

#include "mpibig.h"
#include "taskmessage.h"

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
}  // namespace

void Worker::Run() {
  TaskMessage message;
  do {
    MPI_Status status;
    aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
    MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE, kMainNode,
             kTag, MPI_COMM_WORLD, &status);
    aocommon::SerialIStream stream(std::move(buffer));
    message.Unserialize(stream);

    switch (message.type) {
      case TaskMessage::Type::kGriddingRequest: {
        buffer.resize(message.bodySize);
        MPI_Recv_Big(buffer.data(), message.bodySize, kMainNode, kTag,
                     MPI_COMM_WORLD, &status);
        aocommon::SerialIStream stream(std::move(buffer));
        stream.UInt64();  // skip the nr of packages

        GriddingTask task;
        task.Unserialize(stream);
        scheduler_.Run(std::move(task), [](GriddingResult&) {});
        break;
      }
      case TaskMessage::Type::kLockGrant:
        scheduler_.GrantLock(message.lockId);
        break;
      default:
        break;
    }

  } while (message.type != TaskMessage::Type::kFinish);
  aocommon::Logger::Info << "Worker node " << scheduler_.Rank()
                         << " received exit message.\n";
}
