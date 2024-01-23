#include "worker.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/logger.h>

#include "mpibig.h"
#include "taskmessage.h"

#include "../scheduling/griddingtask.h"
#include "../scheduling/mpischeduler.h"

#include <mpi.h>

#include <cassert>

Worker::Worker(const Settings& settings) : settings_(settings), rank_(-1) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
}

void Worker::Run() {
  TaskMessage message;
  do {
    MPI_Status status;
    aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
    MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE, 0, 0,
             MPI_COMM_WORLD, &status);
    aocommon::SerialIStream stream(std::move(buffer));
    message.Unserialize(stream);

    switch (message.type) {
      case TaskMessage::Type::kGriddingRequest:
        grid(message.bodySize);
        break;
      default:
        break;
    }

  } while (message.type != TaskMessage::Type::kFinish);
  aocommon::Logger::Info << "Worker node " << rank_
                         << " received exit message.\n";
}

void Worker::grid(size_t bodySize) {
  MPI_Status status;
  aocommon::UVector<unsigned char> buffer(bodySize);
  MPI_Recv_Big(buffer.data(), bodySize, 0, 0, MPI_COMM_WORLD, &status);
  aocommon::SerialIStream stream(std::move(buffer));
  stream.UInt64();  // skip the nr of packages

  GriddingTask task;
  task.Unserialize(stream);
  MPIScheduler scheduler(settings_);
  aocommon::Logger::Info << "Worker node " << rank_
                         << " is starting gridding.\n";
  scheduler.RunLocal(std::move(task), [this](GriddingResult& result) {
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

    MPI_Send(msgStream.data(), msgStream.size(), MPI_BYTE, 0, 0,
             MPI_COMM_WORLD);
    MPI_Send_Big(resStream.data(), resStream.size(), 0, 0, MPI_COMM_WORLD);
  });
}
