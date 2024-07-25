#ifndef DIRECT_MS_GRIDDER_H
#define DIRECT_MS_GRIDDER_H

#include <aocommon/image.h>
#include <aocommon/lane.h>

#include "msgridderbase.h"

#include "../structures/resources.h"

namespace wsclean {

class ProgressBar;

template <typename num_t>
class DirectMSGridder final : public MSGridderBase {
 public:
  DirectMSGridder(const Settings& settings, const Resources& resources,
                  MsProviderCollection& ms_provider_collection,
                  size_t gridder_index);

  virtual void Invert() override;

  virtual void Predict(std::vector<aocommon::Image>&& images) override;

  virtual std::vector<aocommon::Image> ResultImages() override {
    return {std::move(_image)};
  }
  virtual size_t getSuggestedWGridSize() const override { return 1; }

 private:
  struct InversionSample {
    num_t uInLambda, vInLambda, wInLambda;
    std::complex<float> sample;
  };
  const Resources _resources;
  aocommon::Image _image;
  aocommon::ImageBase<num_t> _sqrtLMTable;
  std::vector<aocommon::ImageBase<num_t>> _layers;
  aocommon::Lane<InversionSample> _inversionLane;

  void invertMeasurementSet(const MsProviderCollection::MsData& msData,
                            ProgressBar& progress, size_t msIndex);
  void gridSample(const InversionSample& sample, size_t layer);
  aocommon::ImageBase<num_t> GetSqrtLMLookupTable() const;
};

}  // namespace wsclean

#endif
