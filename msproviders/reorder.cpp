#include "reorder.h"

#include <boost/filesystem/path.hpp>

namespace wsclean {
namespace reorder {

std::string GetFilenamePrefix(const std::string& ms_path_str,
                              const std::string& temp_dir) {
  boost::filesystem::path prefix_path;
  if (temp_dir.empty())
    prefix_path = ms_path_str;
  else {
    std::string ms_path_copy(ms_path_str);
    while (!ms_path_copy.empty() && *ms_path_copy.rbegin() == '/')
      ms_path_copy.resize(ms_path_copy.size() - 1);
    boost::filesystem::path ms_path(ms_path_copy);
    prefix_path = boost::filesystem::path(temp_dir) / ms_path.filename();
  }
  std::string prefix(prefix_path.string());
  while (!prefix.empty() && *prefix.rbegin() == '/')
    prefix.resize(prefix.size() - 1);
  return prefix;
}

std::string GetPartPrefix(const std::string& ms_path_str, size_t part_index,
                          aocommon::PolarizationEnum pol, size_t data_desc_id,
                          const std::string& temp_dir) {
  std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream part_prefix;
  part_prefix << prefix << "-part";
  if (part_index < 1000) part_prefix << '0';
  if (part_index < 100) part_prefix << '0';
  if (part_index < 10) part_prefix << '0';
  part_prefix << part_index;
  part_prefix << "-";
  part_prefix << aocommon::Polarization::TypeToShortString(pol);
  part_prefix << "-b" << data_desc_id;
  return part_prefix.str();
}

std::string GetMetaFilename(const std::string& ms_path_str,
                            const std::string& temp_dir, size_t data_desc_id) {
  std::string prefix = GetFilenamePrefix(ms_path_str, temp_dir);

  std::ostringstream s;
  s << prefix << "-spw" << data_desc_id << "-parted-meta.tmp";
  return s.str();
}

size_t GetMaxChannels(const std::vector<ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

std::map<size_t, size_t> GetDataDescIdMap(
    const std::vector<ChannelRange>& channels) {
  std::map<size_t, size_t> data_desc_ids;
  size_t spw_index = 0;
  for (const ChannelRange& range : channels) {
    if (data_desc_ids.count(range.data_desc_id) == 0) {
      data_desc_ids.emplace(range.data_desc_id, spw_index);
      ++spw_index;
    }
  }
  return data_desc_ids;
}

}  // namespace reorder
}  // namespace wsclean
