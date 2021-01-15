#ifndef BUFFERED_MS_GRIDDER_H
#define BUFFERED_MS_GRIDDER_H

#include "../gridding/msgridderbase.h"

#include <memory>

class WGriddingMSGridder : public MSGridderBase {
 public:
  WGriddingMSGridder(size_t threadCount, double memFraction, double absMemLimit,
                     double accuracy);

  virtual void Invert() final override;

  virtual void Predict(ImageF image) final override;
  virtual void Predict(ImageF, ImageF) final override {
    throw std::runtime_error("Can not do imaginary imaging in this mode");
  }

  virtual ImageF ImageRealResult() final override { return std::move(_image); }
  virtual ImageF ImageImaginaryResult() final override {
    throw std::runtime_error("Can not do imaginary imaging in this mode");
  }

  virtual size_t ActualInversionWidth() const final override {
    return _actualInversionWidth;
  }
  virtual size_t ActualInversionHeight() const final override {
    return _actualInversionHeight;
  }

  virtual void FreeImagingData() final override {}

  virtual size_t getSuggestedWGridSize() const final override { return 1; }

 private:
  ImageF _image;

  void gridMeasurementSet(MSData& msData);

  void predictMeasurementSet(MSData& msData);

  size_t calculateMaxNRowsInMemory(size_t channelCount) const;

  void getTrimmedSize(size_t& trimmedWidth, size_t& trimmedHeight) const;

  size_t _cpuCount;
  int64_t _memSize;
  double _accuracy;
  std::unique_ptr<class WGriddingGridder_Simple> _gridder;
};

#endif
