#include "reorderedfilewriter.h"

#include <functional>
#include <stdexcept>

#include <aocommon/logger.h>
#include <aocommon/polarization.h>

using aocommon::Logger;

namespace wsclean {
namespace reordering {

ReorderedFileWriter::ReorderedFileWriter(
    const ReorderedHandleData& data,
    const std::map<size_t, std::set<aocommon::PolarizationEnum>>&
        ms_polarizations_per_data_desc_id,
    double start_time)
    : data_(data),
      ms_polarizations_per_data_desc_id_(ms_polarizations_per_data_desc_id),
      start_time_(start_time),
      channel_parts_(data_.channels_.size()),
      selected_data_desc_ids_(GetDataDescIdMap(data_.channels_)),
      files_(channel_parts_ * data_.polarizations_.size()),
      meta_files_(selected_data_desc_ids_.size()),
      max_channels_(GetMaxChannels(data_.channels_)),
      selected_row_count_per_spw_index_(selected_data_desc_ids_.size(), 0) {
  if (channel_parts_ != 1) {
    Logger::Debug << "Reordering in " << data_.channels_.size() << " channels:";
    for (size_t i = 0; i != channel_parts_; ++i)
      Logger::Debug << ' ' << data_.channels_[i].data_desc_id << ':'
                    << data_.channels_[i].start << '-'
                    << data_.channels_[i].end;
  }
  Logger::Debug << '\n';

  // Each data desc id needs a separate meta file because they can have
  // different uvws and other info.
  size_t file_index = 0;
  for (size_t part = 0; part != channel_parts_; ++part) {
    for (aocommon::PolarizationEnum p : data_.polarizations_) {
      ReorderedDataFiles& file = files_[file_index];
      const std::string part_prefix = GetPartPrefix(
          data_.ms_path_, part, p, data_.channels_[part].data_desc_id,
          data_.temporary_directory_);
      file.data = std::make_unique<std::ofstream>(part_prefix + ".tmp");
      file.weight = std::make_unique<std::ofstream>(part_prefix + "-w.tmp");
      if (data_.initial_model_required_)
        file.model = std::make_unique<std::ofstream>(part_prefix + "-m.tmp");
      file.data->seekp(PartHeader::BINARY_SIZE, std::ios::beg);

      ++file_index;
    }
  }

  // Write header of meta file, one meta file for each data desc id
  // TODO rather than writing we can just skip and write later
  for (const std::pair<const size_t, size_t>& p : selected_data_desc_ids_) {
    const size_t data_desc_id = p.first;
    const size_t spw_index = p.second;
    const std::string meta_filename = GetMetaFilename(
        data_.ms_path_, data_.temporary_directory_, data_desc_id);

    meta_files_[spw_index] = std::make_unique<std::ofstream>(meta_filename);

    reordering::MetaHeader meta_header;
    meta_header.start_time = start_time_;
    meta_header.selected_row_count = 0;  // not yet known
    meta_header.filename_length = data_.ms_path_.size();
    meta_header.Write(*meta_files_[spw_index]);
    meta_files_[spw_index]->write(data_.ms_path_.c_str(),
                                  data_.ms_path_.size());
    if (!meta_files_[spw_index]->good())
      throw std::runtime_error("Error writing to temporary file " +
                               meta_filename);
  }

  selected_rows_total_ = 0;
}

void ReorderedFileWriter::WriteMetaRow(double u, double v, double w,
                                       double time, uint32_t data_desc_id,
                                       uint32_t antenna1, uint32_t antenna2,
                                       uint32_t field_id) {
  reordering::MetaRecord meta;
  meta.u = u;
  meta.v = v;
  meta.w = w;
  meta.time = time;
  meta.antenna1 = antenna1;
  meta.antenna2 = antenna2;
  meta.field_id = field_id;
  const size_t spw_index = selected_data_desc_ids_.find(data_desc_id)->second;
  ++selected_row_count_per_spw_index_[spw_index];
  ++selected_rows_total_;
  std::ofstream& meta_file = *meta_files_[spw_index];
  meta.Write(meta_file);
  if (!meta_file.good())
    throw std::runtime_error("Error writing to temporary file");
}

void ReorderedFileWriter::WriteDataRow(const std::complex<float>* data_array,
                                       const std::complex<float>* model_array,
                                       const float* weight_spectrum_array,
                                       const bool* flag_array,
                                       size_t data_desc_id) {
  const size_t polarizations_per_file =
      aocommon::Polarization::GetVisibilityCount(*data_.polarizations_.begin());
  std::vector<std::complex<float>> data_buffer(polarizations_per_file *
                                               max_channels_);
  std::vector<float> weight_buffer(polarizations_per_file * max_channels_);

  size_t file_index = 0;
  for (size_t part = 0; part != channel_parts_; ++part) {
    if (data_.channels_[part].data_desc_id == data_desc_id) {
      const size_t part_start_ch = data_.channels_[part].start;
      const size_t part_end_ch = data_.channels_[part].end;
      const std::set<aocommon::PolarizationEnum>& ms_polarizations =
          ms_polarizations_per_data_desc_id_.find(data_desc_id)->second;

      for (aocommon::PolarizationEnum p : data_.polarizations_) {
        ReorderedDataFiles& f = files_[file_index];
        reordering::ExtractData(data_buffer.data(), part_start_ch, part_end_ch,
                                ms_polarizations, data_array, p);
        f.data->write(reinterpret_cast<char*>(data_buffer.data()),
                      (part_end_ch - part_start_ch) *
                          sizeof(std::complex<float>) * polarizations_per_file);
        if (!f.data->good())
          throw std::runtime_error("Error writing to temporary data file");

        if (data_.initial_model_required_) {
          reordering::ExtractData(data_buffer.data(), part_start_ch,
                                  part_end_ch, ms_polarizations, model_array,
                                  p);
          f.model->write(reinterpret_cast<char*>(data_buffer.data()),
                         (part_end_ch - part_start_ch) *
                             sizeof(std::complex<float>) *
                             polarizations_per_file);
          if (!f.model->good())
            throw std::runtime_error(
                "Error writing to temporary model data file");
        }

        reordering::ExtractWeights(weight_buffer.data(), part_start_ch,
                                   part_end_ch, ms_polarizations, data_array,
                                   weight_spectrum_array, flag_array, p);
        f.weight->write(reinterpret_cast<char*>(weight_buffer.data()),
                        (part_end_ch - part_start_ch) * sizeof(float) *
                            polarizations_per_file);
        if (!f.weight->good())
          throw std::runtime_error("Error writing to temporary weights file");
        ++file_index;
      }
    } else {
      file_index += data_.polarizations_.size();
    }
  }
}

void ReorderedFileWriter::UpdateMetaHeaders() {
  // Rewrite meta headers to include selected row count
  for (const std::pair<const size_t, size_t>& p : selected_data_desc_ids_) {
    const size_t spw_index = p.second;
    // Data narrowed for reordered file
    reordering::MetaHeader meta_header;
    meta_header.start_time = start_time_;
    meta_header.selected_row_count =
        selected_row_count_per_spw_index_[spw_index];
    meta_header.filename_length = data_.ms_path_.size();
    meta_files_[spw_index]->seekp(0);
    meta_header.Write(*meta_files_[spw_index]);
    meta_files_[spw_index]->write(data_.ms_path_.c_str(),
                                  data_.ms_path_.size());
  }
}

void ReorderedFileWriter::UpdatePartHeaders(bool include_model) {
  size_t file_index = 0;
  for (size_t part = 0; part != channel_parts_; ++part) {
    reordering::PartHeader header;
    header.channel_count =
        data_.channels_[part].end - data_.channels_[part].start;
    header.channel_start = data_.channels_[part].start;
    header.data_desc_id = (uint32_t)data_.channels_[part].data_desc_id;
    header.has_model = include_model;

    for ([[maybe_unused]] const aocommon::PolarizationEnum& pol :
         data_.polarizations_) {
      ReorderedDataFiles& file = files_[file_index];
      file.data->seekp(0, std::ios::beg);
      header.Write(*file.data);
      if (!file.data->good())
        throw std::runtime_error("Error writing to temporary data file");
      ++file_index;
    }
  }
}

void ReorderedFileWriter::PopulateModelWithZeros(
    std::function<void(size_t progress, size_t total)> update_progress) {
  const size_t polarizations_per_file =
      aocommon::Polarization::GetVisibilityCount(*data_.polarizations_.begin());
  std::vector<std::complex<float>> data_buffer(
      polarizations_per_file * max_channels_, {0.0, 0.0});

  for (size_t part = 0; part != channel_parts_; ++part) {
    size_t data_desc_id = data_.channels_[part].data_desc_id;
    size_t channel_count =
        data_.channels_[part].end - data_.channels_[part].start;

    for (const aocommon::PolarizationEnum& pol : data_.polarizations_) {
      std::string part_prefix = reordering::GetPartPrefix(
          data_.ms_path_, part, pol, data_desc_id, data_.temporary_directory_);
      std::ofstream model_file(part_prefix + "-m.tmp");
      const size_t selected_row_count = selected_row_count_per_spw_index_
          [selected_data_desc_ids_.find(data_desc_id)->second];
      for (size_t i = 0; i != selected_row_count; ++i) {
        model_file.write(reinterpret_cast<char*>(data_buffer.data()),
                         channel_count * sizeof(std::complex<float>) *
                             polarizations_per_file);
        update_progress(part * selected_row_count + i,
                        channel_parts_ * selected_row_count);
      }
    }
  }
}

}  // namespace reordering
}  // namespace wsclean
