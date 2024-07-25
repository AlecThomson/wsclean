#include "reorderedmsreader.h"
#include "../reorderedms.h"

namespace wsclean {

ReorderedMsReader::ReorderedMsReader(ReorderedMs* reordered_ms)
    : MSReader(reordered_ms),
      _currentInputRow(0),
      _readPtrRowOffset(0),
      _metaPtrRowOffset(0),
      _weightPtrRowOffset(0) {
  _metaFile.open(ReorderedMs::getMetaFilename(
                     reordered_ms->_handle._data->_msPath,
                     reordered_ms->_handle._data->_temporaryDirectory,
                     reordered_ms->_partHeader.dataDescId),
                 std::ios::in);
  std::vector<char> msPath(reordered_ms->_metaHeader.filenameLength + 1,
                           char(0));
  // meta and data header were read in ReorderedMs constructor
  _metaFile.seekg(ReorderedMs::MetaHeader::BINARY_SIZE, std::ios::beg);
  _metaFile.read(msPath.data(), reordered_ms->_metaHeader.filenameLength);
  std::string partPrefix = ReorderedMs::getPartPrefix(
      msPath.data(), reordered_ms->_partIndex, reordered_ms->_polarization,
      reordered_ms->_partHeader.dataDescId,
      reordered_ms->_handle._data->_temporaryDirectory);
  _dataFile.open(partPrefix + ".tmp", std::ios::in);
  if (!_dataFile.good())
    throw std::runtime_error("Error opening temporary data file in '" +
                             partPrefix + ".tmp'");
  _dataFile.seekg(ReorderedMs::PartHeader::BINARY_SIZE, std::ios::beg);

  _weightFile.open(partPrefix + "-w.tmp", std::ios::in);
  if (!_weightFile.good())
    throw std::runtime_error("Error opening temporary data weight file '" +
                             partPrefix + "-w.tmp'");
}

bool ReorderedMsReader::CurrentRowAvailable() {
  const ReorderedMs& reordered_ms =
      static_cast<const ReorderedMs&>(*_msProvider);
  return _currentInputRow < reordered_ms._metaHeader.selectedRowCount;
}

void ReorderedMsReader::NextInputRow() {
  const ReorderedMs& reordered_ms =
      static_cast<const ReorderedMs&>(*_msProvider);

  ++_currentInputRow;
  if (_currentInputRow < reordered_ms._metaHeader.selectedRowCount) {
    _readPtrRowOffset += 1;
    _metaPtrRowOffset += 1;
    _weightPtrRowOffset += 1;
  }
}

void ReorderedMsReader::ReadMeta(double& u, double& v, double& w) {
  if (_metaPtrRowOffset != 0)
    _metaFile.seekg(_metaPtrRowOffset * (ReorderedMs::MetaRecord::BINARY_SIZE),
                    std::ios::cur);
  _metaPtrRowOffset = -1;

  ReorderedMs::MetaRecord record;
  record.Read(_metaFile);
  u = record.u;
  v = record.v;
  w = record.w;
}

void ReorderedMsReader::ReadMeta(MSProvider::MetaData& metaData) {
  if (_metaPtrRowOffset != 0)
    _metaFile.seekg(_metaPtrRowOffset * (ReorderedMs::MetaRecord::BINARY_SIZE),
                    std::ios::cur);
  _metaPtrRowOffset = -1;

  ReorderedMs::MetaRecord record;
  record.Read(_metaFile);
  metaData.uInM = record.u;
  metaData.vInM = record.v;
  metaData.wInM = record.w;
  metaData.fieldId = record.fieldId;
  metaData.antenna1 = record.antenna1;
  metaData.antenna2 = record.antenna2;
  metaData.time = record.time;
}

void ReorderedMsReader::ReadData(std::complex<float>* buffer) {
  const ReorderedMs& reordered_ms =
      static_cast<const ReorderedMs&>(*_msProvider);

  const int64_t n_visibilities = reordered_ms._partHeader.channelCount *
                                 reordered_ms._polarizationCountInFile;
  if (_readPtrRowOffset != 0) {
    // Data file position was moved forward already, so seek back by one block
    _dataFile.seekg(
        _readPtrRowOffset * (n_visibilities * sizeof(std::complex<float>)),
        std::ios::cur);
  }
  _readPtrRowOffset = -1;
#ifndef NDEBUG
  const size_t pos =
      size_t(_dataFile.tellg()) - ReorderedMs::PartHeader::BINARY_SIZE;
  if (pos != _currentInputRow * n_visibilities * sizeof(std::complex<float>)) {
    std::ostringstream s;
    s << "Not on right pos: " << pos << " instead of "
      << _currentInputRow * n_visibilities * sizeof(std::complex<float>)
      << " (row " << (pos / (n_visibilities * sizeof(std::complex<float>)))
      << " instead of " << _currentInputRow << ")";
    throw std::runtime_error(s.str());
  }
#endif
  _dataFile.read(reinterpret_cast<char*>(buffer),
                 n_visibilities * sizeof(std::complex<float>));
}

void ReorderedMsReader::ReadModel(std::complex<float>* buffer) {
  const ReorderedMs& reordered_ms = static_cast<ReorderedMs&>(*_msProvider);

#ifndef NDEBUG
  if (!reordered_ms._partHeader.hasModel)
    throw std::runtime_error("Reordered MS initialized without model");
#endif
  const size_t rowLength = reordered_ms._partHeader.channelCount *
                           reordered_ms._polarizationCountInFile *
                           sizeof(std::complex<float>);
  std::copy_n(reordered_ms._modelFile.Data() + rowLength * _currentInputRow,
              rowLength, reinterpret_cast<char*>(buffer));
}

void ReorderedMsReader::ReadWeights(float* buffer) {
  const ReorderedMs& reordered_ms =
      static_cast<const ReorderedMs&>(*_msProvider);

  const int64_t n_visibilities = reordered_ms._partHeader.channelCount *
                                 reordered_ms._polarizationCountInFile;
  if (_weightPtrRowOffset != 0) {
    // jump to the previous block of weights
    _weightFile.seekg(_weightPtrRowOffset * (n_visibilities * sizeof(float)),
                      std::ios::cur);
  }
  _weightFile.read(reinterpret_cast<char*>(buffer),
                   n_visibilities * sizeof(float));
  _weightPtrRowOffset = -1;
}

void ReorderedMsReader::WriteImagingWeights(const float* buffer) {
  const ReorderedMs& reordered_ms =
      static_cast<const ReorderedMs&>(*_msProvider);

  if (_imagingWeightsFile == nullptr) {
    std::string partPrefix = ReorderedMs::getPartPrefix(
        reordered_ms._handle._data->_msPath, reordered_ms._partIndex,
        reordered_ms._polarization, reordered_ms._partHeader.dataDescId,
        reordered_ms._handle._data->_temporaryDirectory);
    _imagingWeightsFile.reset(
        new std::fstream(partPrefix + "-imgw.tmp",
                         std::ios::in | std::ios::out | std::ios::binary));
  }
  const size_t nVis = reordered_ms._partHeader.channelCount *
                      reordered_ms._polarizationCountInFile;
  _imagingWeightBuffer.resize(nVis);
  const size_t chunkSize = nVis * sizeof(float);
  _imagingWeightsFile->seekg(chunkSize * _currentInputRow, std::ios::beg);
  _imagingWeightsFile->read(
      reinterpret_cast<char*>(_imagingWeightBuffer.data()),
      nVis * sizeof(float));
  for (size_t i = 0; i != nVis; ++i) {
    if (std::isfinite(buffer[i])) _imagingWeightBuffer[i] = buffer[i];
  }
  _imagingWeightsFile->seekp(chunkSize * _currentInputRow, std::ios::beg);
  _imagingWeightsFile->write(
      reinterpret_cast<const char*>(_imagingWeightBuffer.data()),
      nVis * sizeof(float));
}

}  // namespace wsclean
