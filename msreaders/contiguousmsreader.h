#ifndef MSREADERS_CONTIGUOUSMSREADER_H
#define MSREADERS_CONTIGUOUSMSREADER_H

#include "msreader.h"

class ContiguousMSReader final : public MSReader {
 public:
  ContiguousMSReader(MSProvider* msProvider)
      : MSReader(msProvider),
        _currentInputRow(0),
        _currentInputTimestep(0),
        _currentInputTime(0.0),
        _currentRowId(0),
        _isDataRead(false),
        _isModelRead(false),
        _isWeightRead(false){};
  virtual ~ContiguousMSReader(){};

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
  size_t _currentInputTimestep;
  double _currentInputTime;
  size_t _currentRowId;

  bool _isDataRead, _isModelRead, _isWeightRead;
  std::unique_ptr<casacore::ArrayColumn<float>> _imagingWeightsColumn;

  void readData();

  void readWeights();

  void readModel();
};

#endif