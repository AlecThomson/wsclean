#ifndef INDEPENDENTREADER_H
#define INDEPENDENTREADER_H

#include "msprovider.h"

/**
 * An IndependentReader provides the capability to read (meta)data
 * independent from the location of the reader position in the MSProvider
 */
class IndependentReader {
 public:
  IndependentReader(MSProvider* msProvider) : _msProvider(msProvider){};
  virtual ~IndependentReader(){};

  // TODO: reading (buffered) meta data at the current output position?
  // (potentially lagging behind the _currentInputPosition)
  void GetCachedMetaData(MSProvider::MetaData& metaData) {
    // TODO: this probably is incorrect for a TimeStepBuffer
    // (which doesn't even have a _currentInput/OutputPosition)...

    if ((_msProvider->_currentInputRow - _msProvider->_currentOutputRow) == 1) {
      // Clear buffer as soon as output position catches up with input position?
      _buffer.clear();
      _msProvider->ReadMeta(metaData);
    } else {
      metaData = _buffer[_msProvider->_currentInputRow -
                         _msProvider->_currentOutputRow];
    }
  };
  // TODO: Maybe not even needed, we only need
  // the MetaData
  void GetCachedData(){};

  // BufferMetaData is going to be called during read operations
  void BufferMetaData(const MSProvider::MetaData metaData) {
    if (_msProvider->_currentOutputRow - _msProvider->_currentInputRow <= 1) {
      // Clear buffer, i.e.
      _buffer.clear();
    } else {
      // If _currentInputRow is running ahead of _currentOutputRow, push back to
      // buffer
      _buffer.push_back(metaData);
    }
  };

 private:
  struct RowData {
    // TODO: maybe we don't want to buffer the (model) data
    // but the meta data only
    // std::vector<std::complex<float>> data, model;
    // std::vector<float> weights;
    MSProvider::MetaData metaData;
    size_t rowId;
  };

  MSProvider* _msProvider;

  size_t _bufferPosition;
  std::vector<MSProvider::MetaData> _buffer;
};

#endif
