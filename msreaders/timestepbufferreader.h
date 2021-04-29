#ifndef MSREADERS_TIMESTEPBUFFERREADER_H
#define MSREADERS_TIMESTEPBUFFERREADER_H

#include "msreader.h"
#include "../msproviders/timestepbuffer.h"

class TimestepBufferReader final : public MSReader {
 public:
  TimestepBufferReader(MSProvider* msProvider)
      : MSReader(msProvider), _msReader(msProvider->GetReader().get()) {
    readTimeblock();
  };
  virtual ~TimestepBufferReader(){};

  size_t RowId() const final override { return _buffer[_bufferPosition].rowId; }

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
  void readTimeblock();

  MSReader* _msReader;

  size_t _bufferPosition;
  std::vector<TimestepBuffer::RowData> _buffer;
};

#endif