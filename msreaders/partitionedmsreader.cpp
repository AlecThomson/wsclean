#include "partitionedmsreader.h"
#include "../msproviders/partitionedms.h"

bool PartitionedMSReader::CurrentRowAvailable() {
  const PartitionedMS& partitionedms =
      static_cast<const PartitionedMS&>(*_msProvider);

  return _currentInputRow < partitionedms._metaHeader.selectedRowCount;
}

void PartitionedMSReader::NextInputRow() {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  ++_currentInputRow;
  if (_currentInputRow < partitionedms._metaHeader.selectedRowCount) {
    if (_readPtrIsOk)
      partitionedms._dataFile.seekg(partitionedms._partHeader.channelCount *
                                        partitionedms._polarizationCountInFile *
                                        sizeof(std::complex<float>),
                                    std::ios::cur);
    else
      _readPtrIsOk = true;

    if (_metaPtrIsOk)
      partitionedms._metaFile.seekg(PartitionedMS::MetaRecord::BINARY_SIZE,
                                    std::ios::cur);
    else
      partitionedms._metaPtrIsOk = true;

    if (partitionedms._weightPtrIsOk)
      partitionedms._weightFile.seekg(
          partitionedms._partHeader.channelCount *
              partitionedms._polarizationCountInFile * sizeof(float),
          std::ios::cur);
    _weightPtrIsOk = true;
  }
}

void PartitionedMSReader::ReadMeta(double& u, double& v, double& w,
                                   size_t& dataDescId) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  if (!_metaPtrIsOk)
    partitionedms._metaFile.seekg(-PartitionedMS::MetaRecord::BINARY_SIZE,
                                  std::ios::cur);
  _metaPtrIsOk = false;

  PartitionedMS::MetaRecord record;
  record.read(partitionedms._metaFile);
  u = record.u;
  v = record.v;
  w = record.w;
  dataDescId = record.dataDescId;
}

void PartitionedMSReader::ReadMeta(MSProvider::MetaData& metaData) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  if (!_metaPtrIsOk)
    partitionedms._metaFile.seekg(-PartitionedMS::MetaRecord::BINARY_SIZE,
                                  std::ios::cur);
  _metaPtrIsOk = false;

  PartitionedMS::MetaRecord record;
  record.read(partitionedms._metaFile);
  metaData.uInM = record.u;
  metaData.vInM = record.v;
  metaData.wInM = record.w;
  metaData.dataDescId = record.dataDescId;
  metaData.fieldId = record.fieldId;
  metaData.antenna1 = record.antenna1;
  metaData.antenna2 = record.antenna2;
  metaData.time = record.time;
}

void PartitionedMSReader::ReadData(std::complex<float>* buffer) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  if (!_readPtrIsOk) {
    partitionedms._dataFile.seekg(-partitionedms._partHeader.channelCount *
                                      partitionedms._polarizationCountInFile *
                                      sizeof(std::complex<float>),
                                  std::ios::cur);
  }
#ifndef NDEBUG
  size_t pos = size_t(partitionedms._dataFile.tellg()) -
               sizeof(PartitionedMS::PartHeader);
  if (pos != _currentInputRow * partitionedms._partHeader.channelCount *
                 partitionedms._polarizationCountInFile *
                 sizeof(std::complex<float>)) {
    std::ostringstream s;
    s << "Not on right pos: " << pos << " instead of "
      << _currentInputRow * partitionedms._partHeader.channelCount *
             partitionedms._polarizationCountInFile *
             sizeof(std::complex<float>)
      << " (row "
      << (pos / (partitionedms._partHeader.channelCount *
                 partitionedms._polarizationCountInFile *
                 sizeof(std::complex<float>)))
      << " instead of " << _currentInputRow << ")";
    throw std::runtime_error(s.str());
  }
#endif
  partitionedms._dataFile.read(reinterpret_cast<char*>(buffer),
                               partitionedms._partHeader.channelCount *
                                   partitionedms._polarizationCountInFile *
                                   sizeof(std::complex<float>));
  _readPtrIsOk = false;
}

void PartitionedMSReader::ReadModel(std::complex<float>* buffer) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

#ifndef NDEBUG
  if (!partitionedms._partHeader.hasModel)
    throw std::runtime_error("Partitioned MS initialized without model");
#endif
  size_t rowLength = partitionedms._partHeader.channelCount *
                     partitionedms._polarizationCountInFile *
                     sizeof(std::complex<float>);
  memcpy(
      reinterpret_cast<char*>(buffer),
      partitionedms._modelFileMap + rowLength * partitionedms._currentInputRow,
      rowLength);
}

void PartitionedMSReader::ReadWeights(std::complex<float>* buffer) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  if (!_weightPtrIsOk)
    partitionedms._weightFile.seekg(
        -partitionedms._partHeader.channelCount * sizeof(float), std::ios::cur);
  float* displacedBuffer = reinterpret_cast<float*>(buffer) +
                           partitionedms._partHeader.channelCount *
                               partitionedms._polarizationCountInFile;
  partitionedms._weightFile.read(reinterpret_cast<char*>(displacedBuffer),
                                 partitionedms._partHeader.channelCount *
                                     partitionedms._polarizationCountInFile *
                                     sizeof(float));
  _weightPtrIsOk = false;
  MSProvider::copyRealToComplex(buffer, displacedBuffer,
                                partitionedms._partHeader.channelCount *
                                    partitionedms._polarizationCountInFile);
}

void PartitionedMSReader::ReadWeights(float* buffer) {
  PartitionedMS& partitionedms = static_cast<PartitionedMS&>(*_msProvider);

  if (!_weightPtrIsOk)
    partitionedms._weightFile.seekg(-partitionedms._partHeader.channelCount *
                                        partitionedms._polarizationCountInFile *
                                        sizeof(float),
                                    std::ios::cur);
  partitionedms._weightFile.read(reinterpret_cast<char*>(buffer),
                                 partitionedms._partHeader.channelCount *
                                     partitionedms._polarizationCountInFile *
                                     sizeof(float));
  _weightPtrIsOk = false;
}