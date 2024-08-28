#ifndef MSPROVIDER_REORDERED_FILE_WRITER_
#define MSPROVIDER_REORDERED_FILE_WRITER_

#include "reordering.h"
#include "reorderedhandle.h"

#include <fstream>
#include <cstddef>
#include <vector>
#include <set>
#include <string>
#include <memory>

#include <aocommon/polarization.h>

namespace wsclean {
namespace reordering {

class ReorderedFileWriter {
 private:
  struct ReorderedDataFiles {
    std::unique_ptr<std::ofstream> data;
    std::unique_ptr<std::ofstream> weight;
    std::unique_ptr<std::ofstream> model;
  };

 public:
  ReorderedFileWriter(
      const ReorderedHandleData& data,
      const std::map<size_t, std::set<aocommon::PolarizationEnum>>&
          ms_polarizations_per_data_desc_id,
      double start_time);

  void WriteMetaRow(double u, double v, double w, double time,
                    uint32_t data_desc_id, uint32_t antenna1, uint32_t antenna2,
                    uint32_t field_id);

  void WriteDataRow(const std::complex<float>* data_array,
                    const std::complex<float>* model_array,
                    const float* weight_spectrum_array, const bool* flag_array,
                    size_t data_desc_id);
  void UpdateMetaHeaders();
  void UpdatePartHeaders(bool include_model);
  void PopulateModelWithZeros(
      std::function<void(size_t progress, size_t total)> update_progress);

  ~ReorderedFileWriter() = default;

 private:
  ReorderedHandleData data_;
  std::map<size_t, std::set<aocommon::PolarizationEnum>>
      ms_polarizations_per_data_desc_id_;
  double start_time_;
  size_t channel_parts_;

  // This maps data_desc_id to spw index.
  std::map<size_t, size_t> selected_data_desc_ids_;

  // Ordered as files[pol x channelpart]
  std::vector<ReorderedDataFiles> files_;
  std::vector<std::unique_ptr<std::ofstream>> meta_files_;
  size_t max_channels_;
  aocommon::UVector<size_t> selected_row_count_per_spw_index_;
  size_t selected_rows_total_;
};

}  // namespace reordering
}  // namespace wsclean

#endif
