#ifndef MS_LIST_ITEM_H_
#define MS_LIST_ITEM_H_

#include <string>

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include "../msproviders/msdatadescription.h"

/**
 * Represent an item in a list of (selected) measurement sets. It associates
 * the measurement set with its original filename index on the command line.
 */
struct MsListItem {
  std::unique_ptr<MSDataDescription> ms_description;
  /**
   * Index of the measurement set in the filename list (@ref
   * Settings::filenames). This is useful to associate the measurement set with
   * other information that might depend on the original measurment set index,
   * such as the h5parm solution files (@ref Settings::facetSolutionFiles).
   */
  size_t ms_index;

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.Object(*ms_description).UInt64(ms_index);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    ms_description = MSDataDescription::Unserialize(stream);
    stream.UInt64(ms_index);
  }
};

#endif
