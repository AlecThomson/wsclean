#ifndef WGRIDDING_MS_GRIDDER_H
#define WGRIDDING_MS_GRIDDER_H

#include "../gridding/msgridderbase.h"
#include "../structures/resources.h"

#include <aocommon/image.h>

#include <memory>

#define GRIDDING_PREFETCH_FACTOR 2

class WGriddingGridderBase;

class WGriddingMSGridder final : public MSGridderBase {
 public:
  WGriddingMSGridder(const Settings& settings, const Resources& resources,
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

  void gridMeasurementSet(MSData& msData);

  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getActualTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  const Resources resources_;
  double accuracy_;
  bool use_tuned_wgridder_;
  std::unique_ptr<WGriddingGridderBase> gridder_;

  struct GriddingData {
    bool lastChunk = false;
    bool doneReading = false;
    bool doneComputing = true;
    size_t nRows;
    aocommon::UVector<double> uvwBuffer;
    aocommon::UVector<std::complex<float>> visBuffer;
  };
  GriddingData griddingData[GRIDDING_PREFETCH_FACTOR];

  void read_fn(MSProvider* msProvider, const aocommon::BandData& selectedBand,
               std::complex<float>* modelBuffer, float* weightBuffer,
               bool* isSelected, size_t maxNRows, const size_t dataSize,
               std::vector<std::basic_string<char>>* antennaNames);
};

#endif
