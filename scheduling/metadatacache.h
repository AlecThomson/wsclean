#ifndef META_DATA_CACHE_H
#define META_DATA_CACHE_H

#include <aocommon/io/serialstreamfwd.h>

#include <memory>
#include <vector>

#include "../idg/averagebeam.h"

struct MetaDataCache {
  struct Entry {
    double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM,
        integrationTime;
  };
  std::vector<Entry> msDataVector;
  std::unique_ptr<class AverageBeam> averageBeam;
  // FIXME: find better name?
  // FIXME: assumes one value per call to (de)gridder?!
  std::complex<float> averageBeamCorrection;
  std::complex<float> averageH5Correction;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
