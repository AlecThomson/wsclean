#ifndef TASK_MESSAGE_H
#define TASK_MESSAGE_H

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

struct TaskMessage {
  enum class Type {
    kInvalid,
    kStart,
    kFinish,
    kGriddingRequest,
    kGriddingResult,
  } type;
  union {
    size_t n_writer_groups;  // For kStart type.
    size_t body_size;        // For kGridding* types.
  };

  TaskMessage() : type(Type::kInvalid), body_size(0) {}
  TaskMessage(Type type_, size_t payload) : type(type_), body_size(payload) {}

  constexpr static size_t kSerializedSize = 12;

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(static_cast<std::uint32_t>(type)).UInt64(body_size);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.UInt32(type).UInt64(body_size);
  }
};

#endif
