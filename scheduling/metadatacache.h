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
  /** @{
   *  These variables are incremented with a comparatively small value for each
   * gridded visibility, hence a long double is used to accomodate sufficient
   * precision. The h5_correction_sum is only used when both beam and h5parm
   * solutions are applied. When only one of h5parm or beam solutions are
   * applied, the sum is stored in correction_sum and h5_correction_sum is
   * unused.
   *
   * The correction_sum is always the full combined correction needed for a
   * facet, i.e. the visibilities only need to be corrected by that (combined)
   * correction_sum. However, for correcting the beam, the beam part is
   * recalculated when doing the final image-based beam correction. In order to
   * do so, the h5 sum must be separately stored.
   */
  long double correction_sum = 0.0;
  long double h5_correction_sum = 0.0;
  /**
   * @}
   */
  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
