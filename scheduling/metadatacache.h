#ifndef META_DATA_CACHE_H
#define META_DATA_CACHE_H

#include <aocommon/io/serialstreamfwd.h>

#include <memory>
#include <vector>

struct MetaDataCache {
  struct Entry {
    double min_w;
    double max_w;
    double max_w_with_flags;
    double max_baseline_uvw;
    double max_baseline_in_m;
    double integration_time;
  };
  std::vector<Entry> msDataVector;
  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
