#include "../../msproviders/reordering.h"

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include <aocommon/polarization.h>

#include <vector>
#include <cstddef>

using aocommon::Polarization;
using wsclean::reordering::ChannelRange;
using wsclean::reordering::GetDataDescIdMap;
using wsclean::reordering::GetFilenamePrefix;
using wsclean::reordering::GetMaxChannels;
using wsclean::reordering::GetMetaFilename;
using wsclean::reordering::GetPartPrefix;

BOOST_AUTO_TEST_SUITE(reorder_utils)

BOOST_AUTO_TEST_CASE(filename_prefix) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(filename_prefix_tmp_dir) {
  const std::string filename_prefix = GetFilenamePrefix("tmp/test.ms", "tmp");
  BOOST_CHECK_EQUAL(filename_prefix, "tmp/test.ms");
}

BOOST_AUTO_TEST_CASE(filenameprefix_remove_trailing_separator) {
  const std::string filename_prefix = GetFilenamePrefix("test.ms/", "");
  BOOST_CHECK_EQUAL(filename_prefix, "test.ms");
}

BOOST_AUTO_TEST_CASE(metafilename) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 0);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_with_tmp) {
  const std::string meta_filename = GetMetaFilename("test.ms", "tmp", 0);
  BOOST_CHECK_EQUAL(meta_filename, "tmp/test.ms-spw0-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(metafilename_ddi) {
  const std::string meta_filename = GetMetaFilename("test.ms", "", 1);
  BOOST_CHECK_EQUAL(meta_filename, "test.ms-spw1-parted-meta.tmp");
}

BOOST_AUTO_TEST_CASE(partprefix) {
  const std::string partprefix =
      GetPartPrefix("test.ms", 0, Polarization::StokesI, 0, "");
  BOOST_CHECK_EQUAL(partprefix, "test.ms-part0000-I-b0");
}

BOOST_AUTO_TEST_CASE(max_channel_range) {
  const std::vector<ChannelRange> channel_ranges{
      {0, 0, 100},
      {0, 0, 50},
      {0, 100, 200},
      {0, 50, 500},
  };
  const size_t actual = GetMaxChannels(channel_ranges);
  BOOST_CHECK_EQUAL(actual, 450);
}

BOOST_AUTO_TEST_CASE(max_channel_range_empty_range) {
  const std::vector<ChannelRange> channel_ranges;
  const size_t actual = GetMaxChannels(channel_ranges);
  BOOST_CHECK_EQUAL(actual, 0);
}

BOOST_AUTO_TEST_CASE(data_desc_id_map) {
  const std::vector<ChannelRange> channel_ranges{
      {2, 50, 500},
      {0, 0, 100},
      {1, 0, 50},
  };
  const std::map<size_t, size_t> expected = {{0, 1}, {1, 2}, {2, 0}};
  const std::map<size_t, size_t> actual_ddi_map =
      GetDataDescIdMap(channel_ranges);

  for (const auto& ddi_entry : actual_ddi_map) {
    BOOST_CHECK_EQUAL(ddi_entry.second, expected.find(ddi_entry.first)->second);
  }
}

BOOST_AUTO_TEST_SUITE_END()
