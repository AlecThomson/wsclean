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
  ~WGriddingMSGridder();

  virtual void Invert() override;

  virtual void Predict(std::vector<aocommon::Image>&& images) override;

  virtual std::vector<aocommon::Image> ResultImages() override {
    return {std::move(_image)};
  }

  virtual void FreeImagingData() override {}

  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  aocommon::Image _image;

  std::unique_ptr<WGriddingGridderBase> MakeGridder(size_t width,
                                                    size_t height) const;

  void gridMeasurementSet(const MsProviderCollection::MsData& msData);

  void predictMeasurementSet(const MsProviderCollection::MsData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  const Resources resources_;
  double accuracy_;
  bool use_tuned_wgridder_;
  std::unique_ptr<WGriddingGridderBase> gridder_;
};

}  // namespace wsclean

#endif
