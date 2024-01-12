#ifndef MSPROVIDERS_MSREADERS_REORDEREDMSREADER_H_
#define MSPROVIDERS_MSREADERS_REORDEREDMSREADER_H_

#include "msreader.h"

#include <aocommon/uvector.h>

#include <fstream>

class ReorderedMs;

class ReorderedMsReader final : public MSReader {
 public:
  ReorderedMsReader(ReorderedMs* reordered_ms);
  virtual ~ReorderedMsReader(){};

  size_t RowId() const override { return _currentInputRow; }

  bool CurrentRowAvailable() override;

  void NextInputRow() override;

  void ReadMeta(double& u, double& v, double& w) override;

  void ReadMeta(MSProvider::MetaData& metaData) override;

  void ReadData(std::complex<float>* buffer) override;

  void ReadModel(std::complex<float>* buffer) override;

  void ReadWeights(float* buffer) override;

  void WriteImagingWeights(const float* buffer) override;

 private:
  size_t _currentInputRow;

  // Chunkoffset counts the amount of data rows we are ahead or behind.
  // Positive values mean we are lagging, whereas negative values mean we are
  //  ahead of the current time step.
  long _readPtrRowOffset, _metaPtrRowOffset, _weightPtrRowOffset;

  std::ifstream _metaFile, _weightFile, _dataFile;

  aocommon::UVector<float> _imagingWeightBuffer;
  std::unique_ptr<std::fstream> _imagingWeightsFile;
};

#endif
