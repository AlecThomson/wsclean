#ifndef MSPROVIDERS_REORDERED_MS_PROVIDER_
#define MSPROVIDERS_REORDERED_MS_PROVIDER_

#include "msprovider.h"
#include "reorder.h"

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

namespace wsclean {

class ReorderedMsReader;

class ReorderedMsProvider final : public MSProvider {
  friend class ReorderedMsReader;

 public:
  class Handle;

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
                          const std::vector<reorder::ChannelRange>& channels,
                          const class MSSelection& selection,
                          const string& data_column_name, bool include_model,
                          bool initial_model_required,
                          const class Settings& settings);

  static Handle GenerateHandleFromReorderedData(
      const std::string& ms_path, const string& data_column_name,
      const std::string& temporary_directory,
      const std::vector<reorder::ChannelRange>& channels,
      bool initial_model_required, bool model_update_required,
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
                 const std::vector<reorder::ChannelRange>& channels,
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
      std::vector<reorder::ChannelRange> channels_;
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
           const std::vector<reorder::ChannelRange>& channels,
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

  const Handle handle_;
  const size_t part_index_;
  const size_t data_desc_id_;
  MappedFile model_file_;
  size_t current_output_row_;
  std::unique_ptr<std::ofstream> model_data_file_;
  const aocommon::PolarizationEnum polarization_;
  size_t polarization_count_in_file_;

  reorder::MetaHeader meta_header_;
  reorder::PartHeader part_header_;
};

}  // namespace wsclean

#endif
