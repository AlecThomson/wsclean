#ifndef INDEPENDENTREADER_H
#define INDEPENDENTREADER_H

#include "msprovider.h"

/**
 * An IndependentReader provides the capability to read (meta)data
 * independent from the location of the reader position in the MSProvider
 */
class IndependentReader {
 public:
  IndependentReader(const MSProvider* msProvider) : _msProvider(msProvider){};
  virtual ~IndependentReader(){};
  // TODO: reading the meta data at the current output position?!
  void ReadMetaAtWritePosition(){};

 private:
  struct RowData {
    std::vector<std::complex<float>> data, model;
    std::vector<float> weights;
    MSProvider::MetaData metaData;
    size_t rowId;
  };

  const MSProvider* _msProvider;

  size_t _bufferPosition;
  std::vector<RowData> _buffer;
};

#endif
