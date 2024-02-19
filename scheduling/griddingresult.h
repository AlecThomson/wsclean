#ifndef GRIDDING_RESULT_H
#define GRIDDING_RESULT_H

#include <aocommon/image.h>
#include <aocommon/io/serialstreamfwd.h>

#include <string>

#include "metadatacache.h"

class AverageBeam;

struct GriddingResult {
  GriddingResult();
  GriddingResult(const GriddingResult& source) = delete;
  GriddingResult(GriddingResult&& source) noexcept;
  ~GriddingResult();
  GriddingResult& operator=(const GriddingResult& rhs) = delete;
  GriddingResult& operator=(GriddingResult&& rhs) noexcept;

  // The unique_id of the GriddingTask that produced the result.
  uint32_t unique_id = 0;
  // These values are equal for all facets.
  double startTime = 0.0;
  // There could be reasons to make the beam size facet specific. In the
  // current implementation, it is the same for all facets.
  double beamSize = 0.0;
  size_t griddedVisibilityCount = 0;
  double visibilityWeightSum = 0.0;

  /// FacetData holds values that differ between facets.
  struct FacetData {
    /**
     * List of produced images. When performing complex images, images[0] will
     * be the real part and images[1] will be the imaginary part. When
     * performing full polarization imaging (indicated with
     * Polarization::FullStokes) with IDG, the list will contain all four images
     * ordered IQUV. In all other cases, this list will only hold one image.
     */
    std::vector<aocommon::Image> images;

    double imageWeight = 0.0;
    double normalizationFactor = 0.0;
    size_t actualWGridSize = 0;
    double effectiveGriddedVisibilityCount = 0.0;
    /// See VisibilityModifier for an explanation of these parameters
    double averageCorrection = 0.0;
    double averageH5Correction = 0.0;

    std::unique_ptr<MetaDataCache> cache;
    std::unique_ptr<AverageBeam> averageBeam;

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);
  };
  std::vector<FacetData> facets;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
