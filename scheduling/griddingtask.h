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
  bool verbose;
  std::unique_ptr<AverageBeam> averageBeam;
  bool storeImagingWeights;

  std::shared_ptr<ImageWeights> imageWeights;
  std::vector<std::unique_ptr<MSDataDescription>> msList;

  /**
   * Images for prediction. See the documentation of
   * @ref GriddingResult::images for an explanation of why this is a vector.
   */
  std::vector<aocommon::Image> modelImages;
  ObservationInfo observationInfo;

  struct FacetData {
    /// The default constructor is required since deserializing an ObjectVector
    /// first creates empty objects and then calls Unserialize on those objects.
    FacetData() = default;

    explicit FacetData(
        size_t _index, double _l_shift, double _m_shift,
        std::unique_ptr<MetaDataCache> _cache,
        const std::shared_ptr<schaapcommon::facets::Facet>& _facet)
        : index(_index),
          l_shift(_l_shift),
          m_shift(_m_shift),
          cache(std::move(_cache)),
          facet(_facet) {}

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    size_t index;    ///< Index of the facet, between zero and n_facets.
    double l_shift;  ///< l_shift, adjusted to the center of the facet.
    double m_shift;  ///< m_shift, adjusted to the center of the facet.
    std::unique_ptr<MetaDataCache> cache;
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
