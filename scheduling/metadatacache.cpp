#include "metadatacache.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void MetaDataCache::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt64(msDataVector.size());
  for (const Entry& entry : msDataVector) {
    stream.Double(entry.min_w)
        .Double(entry.max_w)
        .Double(entry.max_w_with_flags)
        .Double(entry.max_baseline_uvw)
        .Double(entry.max_baseline_in_m)
        .Double(entry.integration_time);
  }
}

void MetaDataCache::Unserialize(aocommon::SerialIStream& stream) {
  msDataVector.resize(stream.UInt64());
  for (Entry& entry : msDataVector) {
    stream.Double(entry.min_w)
        .Double(entry.max_w)
        .Double(entry.max_w_with_flags)
        .Double(entry.max_baseline_uvw)
        .Double(entry.max_baseline_in_m)
        .Double(entry.integration_time);
  }
}
