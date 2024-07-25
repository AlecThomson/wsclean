#include "timestepbuffer.h"
#include "msreaders/timestepbufferreader.h"

namespace wsclean {

std::unique_ptr<MSReader> TimestepBuffer::MakeReader() {
  std::unique_ptr<MSReader> reader(new TimestepBufferReader(this));
  return reader;
}

}  // namespace wsclean
