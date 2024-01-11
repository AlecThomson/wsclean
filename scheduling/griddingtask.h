#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"

#include "../msproviders/msdatadescription.h"

#include "metadatacache.h"

class AverageBeam;

namespace schaapcommon {
namespace facets {
class Facet;
}
}  // namespace schaapcommon

class GriddingTask {
 public:
  GriddingTask();
  GriddingTask(const GriddingTask&) = delete;
  GriddingTask(GriddingTask&& source) noexcept;
  ~GriddingTask() noexcept;
  GriddingTask& operator=(const GriddingTask& source) = delete;
  GriddingTask& operator=(GriddingTask&& source) noexcept;

  enum Operation { Invert, Predict } operation;
  bool imagePSF;
  bool subtractModel;
  aocommon::PolarizationEnum polarization;
  bool isFirstTask;
  bool storeImagingWeights;

  std::shared_ptr<ImageWeights> imageWeights;
  std::vector<std::unique_ptr<MSDataDescription>> msList;

  ObservationInfo observationInfo;

  struct FacetData {
    /// The default constructor is required since deserializing an ObjectVector
    /// first creates empty objects and then calls Unserialize on those objects.
    FacetData() = default;

    explicit FacetData(
        size_t _index, double _l_shift, double _m_shift,
        std::unique_ptr<MetaDataCache> _cache,
        const std::shared_ptr<schaapcommon::facets::Facet>& _facet);

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    /// Images for prediction. The documentation of
    /// @ref GriddingResult::FacetData::images explains why it is a vector.
    std::vector<aocommon::Image> modelImages;

    size_t index;    ///< Index of the facet, between zero and n_facets.
    double l_shift;  ///< l_shift, adjusted to the center of the facet.
    double m_shift;  ///< m_shift, adjusted to the center of the facet.
    std::unique_ptr<MetaDataCache> cache;
    std::unique_ptr<AverageBeam> averageBeam;
    /// The facet itself. If null, faceting is disabled.
    std::shared_ptr<schaapcommon::facets::Facet> facet;
  };

  /// 'facets' always contains at least one element.
  /// When faceting is disabled, there is always a single element with
  /// a null 'facet' pointer. That 'facet' covers the entire image.
  std::vector<FacetData> facets;
  size_t facetGroupIndex;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
