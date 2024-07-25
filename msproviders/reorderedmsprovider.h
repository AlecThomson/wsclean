#ifndef MSPROVIDERS_REORDERED_MS_PROVIDER_
#define MSPROVIDERS_REORDERED_MS_PROVIDER_

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

namespace wsclean {

class ReorderedMsReader;

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
  uint16_t field_id;
};
struct PartHeaderBuffer {
  uint64_t channel_count;
  uint64_t channel_start;
  uint32_t data_desc_id;
  bool has_model;
};
struct MetaHeaderBuffer {
  double start_time;
  uint64_t selected_row_count;
  uint32_t filename_length;
};
#pragma pack(pop)

class ReorderedMsProvider final : public MSProvider {
  friend class ReorderedMsReader;

 public:
  class Handle;

  struct ChannelRange {
    int data_desc_id;
    size_t start, end;
    bool operator<(const ChannelRange& rhs) const {
      if (data_desc_id < rhs.data_desc_id) return true;
      if (data_desc_id > rhs.data_desc_id) return false;
      if (start < rhs.start) return true;
      if (start > rhs.start) return false;
      return end < rhs.end;
    }
  };

  ReorderedMsProvider(const Handle& handle, size_t part_index,
                      aocommon::PolarizationEnum polarization,
                      size_t data_desc_id);

  virtual ~ReorderedMsProvider();

  ReorderedMsProvider(const ReorderedMsProvider&) = delete;
  ReorderedMsProvider& operator=(const ReorderedMsProvider&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  SynchronizedMS MS() override {
    return SynchronizedMS(handle_.data_->ms_path_.data());
  }

  const std::string& DataColumnName() override {
    return handle_.data_->data_column_name_;
  }

  void NextOutputRow() override;

  void ResetWritePosition() override { current_output_row_ = 0; };

  void WriteModel(const std::complex<float>* buffer, bool add_to_MS) override;

  void ReopenRW() override {}

  double StartTime() override { return meta_header_.start_time; }

  void MakeIdToMSRowMapping(std::vector<size_t>& id_to_MS_row) override;

  aocommon::PolarizationEnum Polarization() override { return polarization_; }

  size_t NChannels() override { return part_header_.channel_count; }
  size_t NPolarizations() override { return polarization_count_in_file_; }
  size_t NAntennas() override { return handle_.data_->n_antennas_; }

  size_t DataDescId() override { return part_header_.data_desc_id; }

  static Handle Partition(const string& ms_path,
                          const std::vector<ChannelRange>& channels,
                          const class MSSelection& selection,
                          const string& data_column_name, bool include_model,
                          bool initial_model_required,
                          const class Settings& settings);

  static Handle GenerateHandleFromReorderedData(
      const std::string& ms_path, const string& data_column_name,
      const std::string& temporary_directory,
      const std::vector<ChannelRange>& channels, bool initial_model_required,
      bool model_update_required,
      const std::set<aocommon::PolarizationEnum>& polarizations,
      const MSSelection& selection, const aocommon::MultiBandData& bands,
      size_t n_antennas, bool keep_temporary_files);

  const aocommon::BandData& Band() override {
    return handle_.data_->bands_[data_desc_id_];
  }

  class Handle {
    // ReorderedMsReader is a friend of Handle
    // in order to access the data_ member.
    friend class ReorderedMsReader;

   public:
    Handle() = default;

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

    friend class ReorderedMsProvider;

   private:
    struct HandleData {
      HandleData() : is_copy_(false) {}

      HandleData(const std::string& ms_path, const string& data_column_name,
                 const std::string& temporary_directory,
                 const std::vector<ChannelRange>& channels,
                 bool initial_model_required, bool model_update_required,
                 const std::set<aocommon::PolarizationEnum>& polarizations,
                 const MSSelection& selection,
                 const aocommon::MultiBandData& bands, size_t n_antennas,
                 bool keep_temporary_files)
          : ms_path_(ms_path),
            data_column_name_(data_column_name),
            temporary_directory_(temporary_directory),
            channels_(channels),
            initial_model_required_(initial_model_required),
            model_update_required_(model_update_required),
            polarizations_(polarizations),
            selection_(selection),
            bands_(bands),
            n_antennas_(n_antennas),
            is_copy_(false),
            keep_temporary_files_(keep_temporary_files) {}

      ~HandleData();

      std::string ms_path_, data_column_name_, temporary_directory_;
      std::vector<ChannelRange> channels_;
      bool initial_model_required_, model_update_required_;
      std::set<aocommon::PolarizationEnum> polarizations_;
      MSSelection selection_;
      aocommon::MultiBandData bands_;
      size_t n_antennas_;
      bool is_copy_;
      bool keep_temporary_files_;

      void Serialize(aocommon::SerialOStream& stream) const;
      void Unserialize(aocommon::SerialIStream& stream);
    };
    std::shared_ptr<HandleData> data_;

    Handle(const std::string& ms_path, const string& data_column_name,
           const std::string& temporary_directory,
           const std::vector<ChannelRange>& channels,
           bool initial_model_required, bool model_update_required,
           const std::set<aocommon::PolarizationEnum>& polarizations,
           const MSSelection& selection, const aocommon::MultiBandData& bands,
           size_t n_antennas, bool keep_temporary_files)
        : data_(std::make_shared<HandleData>(
              ms_path, data_column_name, temporary_directory, channels,
              initial_model_required, model_update_required, polarizations,
              selection, bands, n_antennas, keep_temporary_files)) {}
  };

 private:
  static void unpartition(const Handle::HandleData& handle);

  /**
   * Make a map that maps data_desc_id to spw (spectral window) index.
   */
  static std::map<size_t, size_t> GetDataDescIdMap(
      const std::vector<ChannelRange>& channels);

  const Handle handle_;
  const size_t part_index_;
  const size_t data_desc_id_;
  MappedFile model_file_;
  size_t current_output_row_;
  std::unique_ptr<std::ofstream> model_data_file_;
  const aocommon::PolarizationEnum polarization_;
  size_t polarization_count_in_file_;

  struct MetaHeader {
    double start_time = 0.0;
    uint64_t selected_row_count = 0;
    uint32_t filename_length = 0;
    void Read(std::istream& str) {
      MetaHeaderBuffer meta_header_buffer;
      str.read(reinterpret_cast<char*>(&meta_header_buffer),
               sizeof(MetaHeaderBuffer));
      start_time = meta_header_buffer.start_time;
      selected_row_count = meta_header_buffer.selected_row_count;
      filename_length = meta_header_buffer.filename_length;
    }
    void Write(std::ostream& str) const {
      MetaHeaderBuffer meta_header_buffer{start_time, selected_row_count,
                                          filename_length};
      str.write(reinterpret_cast<const char*>(&meta_header_buffer),
                sizeof(MetaHeaderBuffer));
    }
    static constexpr size_t BINARY_SIZE = sizeof(start_time) +
                                          sizeof(selected_row_count) +
                                          sizeof(filename_length);
    static_assert(BINARY_SIZE == 20);
  } meta_header_;
  struct MetaRecord {
    double u = 0.0, v = 0.0, w = 0.0, time = 0.0;
    uint16_t antenna1 = 0, antenna2 = 0, field_id = 0;
    static constexpr size_t BINARY_SIZE =
        sizeof(double) * 4 + sizeof(uint16_t) * 3;
    static_assert(BINARY_SIZE == 38);
    void Read(std::istream& str) {
      MetaRecordBuffer meta_record_buffer;
      str.read(reinterpret_cast<char*>(&meta_record_buffer),
               sizeof(MetaRecordBuffer));

      u = meta_record_buffer.u;
      v = meta_record_buffer.v;
      w = meta_record_buffer.w;
      time = meta_record_buffer.time;
      antenna1 = meta_record_buffer.antenna1;
      antenna2 = meta_record_buffer.antenna2;
      field_id = meta_record_buffer.field_id;
    }
    void Write(std::ostream& str) const {
      MetaRecordBuffer meta_record_buffer{u,        v,        w,       time,
                                          antenna1, antenna2, field_id};
      str.write(reinterpret_cast<const char*>(&meta_record_buffer),
                sizeof(MetaRecordBuffer));
    }
  };
  struct PartHeader {
    uint64_t channel_count = 0;
    uint64_t channel_start = 0;
    uint32_t data_desc_id = 0;
    bool has_model = false;
    static constexpr size_t BINARY_SIZE =
        sizeof(channel_count) + sizeof(channel_start) + sizeof(data_desc_id) +
        sizeof(has_model);
    static_assert(BINARY_SIZE == 21);
    void Read(std::istream& str) {
      PartHeaderBuffer part_header_buffer;
      str.read(reinterpret_cast<char*>(&part_header_buffer),
               sizeof(PartHeaderBuffer));
      channel_count = part_header_buffer.channel_count;
      channel_start = part_header_buffer.channel_start;
      data_desc_id = part_header_buffer.data_desc_id;
      has_model = part_header_buffer.has_model;
    }
    void Write(std::ostream& str) const {
      PartHeaderBuffer part_header_buffer{channel_count, channel_start,
                                          data_desc_id, has_model};
      str.write(reinterpret_cast<const char*>(&part_header_buffer),
                sizeof(PartHeaderBuffer));
    }
  } part_header_;

  static std::string GetFilenamePrefix(const std::string& ms_path,
                                       const std::string& temp_dir);
  static std::string GetPartPrefix(const std::string& ms_path,
                                   size_t part_index,
                                   aocommon::PolarizationEnum pol,
                                   size_t data_desc_id,
                                   const std::string& temp_dir);
  static std::string GetMetaFilename(const std::string& ms_path,
                                     const std::string& temp_dir,
                                     size_t data_desc_id);
};

}  // namespace wsclean

#endif
