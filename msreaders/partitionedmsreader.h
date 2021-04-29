#ifndef MSREADERS_PARTITIONEDMSREADER_H
#define MSREADERS_PARTITIONEDMSREADER_H

#include "msreader.h"

class PartitionedMSReader final : public MSReader {
 public:
  PartitionedMSReader(MSProvider* msProvider)
      : MSReader(msProvider),
        _currentInputRow(0),
        _readPtrIsOk(true),
        _metaPtrIsOk(true),
        _weightPtrIsOk(true){};
  virtual ~PartitionedMSReader(){};

  bool CurrentRowAvailable() final override;

  void NextInputRow() final override;

  void ReadMeta(double& u, double& v, double& w,
                size_t& dataDescId) final override;

  void ReadMeta(MSProvider::MetaData& metaData) final override;

  void ReadData(std::complex<float>* buffer) final override;

  void ReadModel(std::complex<float>* buffer) final override;

  void ReadWeights(std::complex<float>* buffer) final override;

  void ReadWeights(float* buffer) final override;

  void WriteImagingWeights(const float* buffer) final override;

 private:
  size_t _currentInputRow;
  bool _readPtrIsOk, _metaPtrIsOk, _weightPtrIsOk;

  // FIXME: why is _modelDataFile needed?
  std::unique_ptr<std::ofstream> _modelDataFile;
  std::unique_ptr<std::fstream> _imagingWeightsFile;
};

#endif