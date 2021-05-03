#ifndef MSREADERS_PARTITIONEDMSREADER_H
#define MSREADERS_PARTITIONEDMSREADER_H

#include "msreader.h"

#include <aocommon/uvector.h>

#include <fstream>

class PartitionedMS;

class PartitionedMSReader final : public MSReader {
 public:
  PartitionedMSReader(PartitionedMS* partitionedMS);
  virtual ~PartitionedMSReader(){};

  size_t RowId() const final override { return _currentInputRow; }

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

  std::ifstream _metaFile, _weightFile, _dataFile;

  // FIXME: both _weightBuffer and _modelBuffer seem redundant
  // and have been removed. Check with Andre
  aocommon::UVector<float> _imagingWeightBuffer;

  // FIXME: _modelDataFile seems redundant (here and in
  // PartitionedMS class)
  std::unique_ptr<std::ofstream> _modelDataFile;
  std::unique_ptr<std::fstream> _imagingWeightsFile;
};

#endif