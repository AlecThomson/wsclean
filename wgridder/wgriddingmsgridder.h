#ifndef WGRIDDING_MS_GRIDDER_H
#define WGRIDDING_MS_GRIDDER_H

#include "../gridding/msgridderbase.h"
#include "../structures/resources.h"

#include <aocommon/image.h>

#include <memory>

namespace wsclean {

class WGriddingGridderBase;

class WGriddingMSGridder final : public MSGridderBase {
 public:
  WGriddingMSGridder(const Settings& settings, const Resources& resources,
                     MsProviderCollection& ms_provider_collection,
                     bool use_tuned_wgridder);
  ~WGriddingMSGridder() final;

  void StartInversion() final;
  size_t GridMeasurementSet(const MsProviderCollection::MsData& ms_data) final;
  void FinishInversion() final;

  void StartPredict(std::vector<aocommon::Image>&& images) final;
  size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) final;
  void FinishPredict() final;

  std::vector<aocommon::Image> ResultImages() final {
    return {std::move(_image)};
  }

  void FreeImagingData() final {}

  size_t getSuggestedWGridSize() const final { return 1; }

 private:
  aocommon::Image _image;

  std::unique_ptr<WGriddingGridderBase> MakeGridder(size_t width,
                                                    size_t height) const;

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  const Resources resources_;
  double accuracy_;
  bool use_tuned_wgridder_;
  std::unique_ptr<WGriddingGridderBase> gridder_;
};

}  // namespace wsclean

#endif
