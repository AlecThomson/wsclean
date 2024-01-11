#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include "msprovider.h"

#include "../structures/msselection.h"
#include "../system/mappedfile.h"

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <fstream>
#include <string>
#include <map>

class PartitionedMSReader;

// We will create some efficiently packed structs to fetch data with 1 read.
// This will reduce the count of file-reads that are made.
// We do not apply this on the classes/structs themselves, because this will
//  reduce the data access performance.
#pragma pack(push, 1)
struct MetaRecordBuffer {
  double u;
  double v;
  double w;
  double time;
  uint16_t antenna1;
  uint16_t antenna2;
  uint16_t fieldId;
};
struct PartHeaderBuffer {
  uint64_t channelCount;
  uint64_t channelStart;
  uint32_t dataDescId;
  bool hasModel;
};
struct MetaHeaderBuffer {
  double startTime;
  uint64_t selectedRowCount;
  uint32_t filenameLength;
};
#pragma pack(pop)

class PartitionedMS final : public MSProvider {
  friend class PartitionedMSReader;

 public:
  class Handle;

  struct ChannelRange {
    int dataDescId;
    size_t start, end;
    bool operator<(const ChannelRange& rhs) const {
      if (dataDescId < rhs.dataDescId) return true;
      if (dataDescId > rhs.dataDescId) return false;
      if (start < rhs.start) return true;
      if (start > rhs.start) return false;
      return end < rhs.end;
    }
  };

  PartitionedMS(const Handle& handle, size_t partIndex,
                aocommon::PolarizationEnum polarization, size_t dataDescId);

  virtual ~PartitionedMS();

  PartitionedMS(const PartitionedMS&) = delete;
  PartitionedMS& operator=(const PartitionedMS&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  SynchronizedMS MS() override {
    return SynchronizedMS(_handle._data->_msPath.data());
  }

  const std::string& DataColumnName() override {
    return _handle._data->_dataColumnName;
  }

  void NextOutputRow() override;

  void ResetWritePosition() override { _currentOutputRow = 0; };

  void WriteModel(const std::complex<float>* buffer, bool addToMS) override;

  void ReopenRW() override {}

  double StartTime() override { return _metaHeader.startTime; }

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override;

  aocommon::PolarizationEnum Polarization() override { return _polarization; }

  size_t NChannels() override { return _partHeader.channelCount; }
  size_t NPolarizations() override { return _polarizationCountInFile; }
  size_t NAntennas() override { return _handle._data->_nAntennas; }

  size_t DataDescId() override { return _partHeader.dataDescId; }

  static Handle Partition(const string& msPath,
                          const std::vector<ChannelRange>& channels,
                          const class MSSelection& selection,
                          const string& dataColumnName, bool includeModel,
                          bool initialModelRequired,
                          const class Settings& settings);

  const aocommon::BandData& Band() override {
    return _handle._data->_bands[_dataDescId];
  }

  class Handle {
    // PartitionedMSReader is a friend of Handle
    // in order to access the _data member.
    friend class PartitionedMSReader;

   public:
    Handle() = default;

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    friend class PartitionedMS;

   private:
    struct HandleData {
      HandleData() : _isCopy(false) {}

      HandleData(const std::string& msPath, const string& dataColumnName,
                 const std::string& temporaryDirectory,
                 const std::vector<ChannelRange>& channels,
                 bool initialModelRequired, bool modelUpdateRequired,
                 const std::set<aocommon::PolarizationEnum>& polarizations,
                 const MSSelection& selection,
                 const aocommon::MultiBandData& bands, size_t nAntennas)
          : _msPath(msPath),
            _dataColumnName(dataColumnName),
            _temporaryDirectory(temporaryDirectory),
            _channels(channels),
            _initialModelRequired(initialModelRequired),
            _modelUpdateRequired(modelUpdateRequired),
            _polarizations(polarizations),
            _selection(selection),
            _bands(bands),
            _nAntennas(nAntennas),
            _isCopy(false) {}

      ~HandleData();

      std::string _msPath, _dataColumnName, _temporaryDirectory;
      std::vector<ChannelRange> _channels;
      bool _initialModelRequired, _modelUpdateRequired;
      std::set<aocommon::PolarizationEnum> _polarizations;
      MSSelection _selection;
      aocommon::MultiBandData _bands;
      size_t _nAntennas;
      bool _isCopy;

      void Serialize(aocommon::SerialOStream& stream) const;
      void Unserialize(aocommon::SerialIStream& stream);
    };
    std::shared_ptr<HandleData> _data;

    Handle(const std::string& msPath, const string& dataColumnName,
           const std::string& temporaryDirectory,
           const std::vector<ChannelRange>& channels, bool initialModelRequired,
           bool modelUpdateRequired,
           const std::set<aocommon::PolarizationEnum>& polarizations,
           const MSSelection& selection, const aocommon::MultiBandData& bands,
           size_t nAntennas)
        : _data(std::make_shared<HandleData>(
              msPath, dataColumnName, temporaryDirectory, channels,
              initialModelRequired, modelUpdateRequired, polarizations,
              selection, bands, nAntennas)) {}
  };

 private:
  static void unpartition(const Handle::HandleData& handle);

  /**
   * Make a map that maps dataDescId to spw (spectral window) index.
   */
  static std::map<size_t, size_t> getDataDescIdMap(
      const std::vector<ChannelRange>& channels);

  const Handle _handle;
  const size_t _partIndex;
  const size_t _dataDescId;
  MappedFile _modelFile;
  size_t _currentOutputRow;
  std::unique_ptr<std::ofstream> _modelDataFile;
  const aocommon::PolarizationEnum _polarization;
  size_t _polarizationCountInFile;

  struct MetaHeader {
    double startTime = 0.0;
    uint64_t selectedRowCount = 0;
    uint32_t filenameLength = 0;
    void Read(std::istream& str) {
      MetaHeaderBuffer metaHeaderBuffer;
      str.read(reinterpret_cast<char*>(&metaHeaderBuffer),
               sizeof(MetaHeaderBuffer));
      startTime = metaHeaderBuffer.startTime;
      selectedRowCount = metaHeaderBuffer.selectedRowCount;
      filenameLength = metaHeaderBuffer.filenameLength;
    }
    void Write(std::ostream& str) const {
      MetaHeaderBuffer metaHeaderBuffer{startTime, selectedRowCount,
                                        filenameLength};
      str.write(reinterpret_cast<const char*>(&metaHeaderBuffer),
                sizeof(MetaHeaderBuffer));
    }
    static constexpr size_t BINARY_SIZE =
        sizeof(startTime) + sizeof(selectedRowCount) + sizeof(filenameLength);
    static_assert(BINARY_SIZE == 20);
  } _metaHeader;
  struct MetaRecord {
    double u = 0.0, v = 0.0, w = 0.0, time = 0.0;
    uint16_t antenna1 = 0, antenna2 = 0, fieldId = 0;
    static constexpr size_t BINARY_SIZE =
        sizeof(double) * 4 + sizeof(uint16_t) * 3;
    static_assert(BINARY_SIZE == 38);
    void Read(std::istream& str) {
      MetaRecordBuffer metaRecordBuffer;
      str.read(reinterpret_cast<char*>(&metaRecordBuffer),
               sizeof(MetaRecordBuffer));

      u = metaRecordBuffer.u;
      v = metaRecordBuffer.v;
      w = metaRecordBuffer.w;
      time = metaRecordBuffer.time;
      antenna1 = metaRecordBuffer.antenna1;
      antenna2 = metaRecordBuffer.antenna2;
      fieldId = metaRecordBuffer.fieldId;
    }
    void Write(std::ostream& str) const {
      MetaRecordBuffer metaRecordBuffer{u,        v,        w,      time,
                                        antenna1, antenna2, fieldId};
      str.write(reinterpret_cast<const char*>(&metaRecordBuffer),
                sizeof(MetaRecordBuffer));
    }
  };
  struct PartHeader {
    uint64_t channelCount = 0;
    uint64_t channelStart = 0;
    uint32_t dataDescId = 0;
    bool hasModel = false;
    static constexpr size_t BINARY_SIZE = sizeof(channelCount) +
                                          sizeof(channelStart) +
                                          sizeof(dataDescId) + sizeof(hasModel);
    static_assert(BINARY_SIZE == 21);
    void Read(std::istream& str) {
      PartHeaderBuffer partHeaderBuffer;
      str.read(reinterpret_cast<char*>(&partHeaderBuffer),
               sizeof(PartHeaderBuffer));
      channelCount = partHeaderBuffer.channelCount;
      channelStart = partHeaderBuffer.channelStart;
      dataDescId = partHeaderBuffer.dataDescId;
      hasModel = partHeaderBuffer.hasModel;
    }
    void Write(std::ostream& str) const {
      PartHeaderBuffer partHeaderBuffer{channelCount, channelStart, dataDescId,
                                        hasModel};
      str.write(reinterpret_cast<const char*>(&partHeaderBuffer),
                sizeof(PartHeaderBuffer));
    }
  } _partHeader;

  static std::string getFilenamePrefix(const std::string& msPath,
                                       const std::string& tempDir);
  static std::string getPartPrefix(const std::string& msPath, size_t partIndex,
                                   aocommon::PolarizationEnum pol,
                                   size_t dataDescId,
                                   const std::string& tempDir);
  static std::string getMetaFilename(const std::string& msPath,
                                     const std::string& tempDir,
                                     size_t dataDescId);
};

#endif
