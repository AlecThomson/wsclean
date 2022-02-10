#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  entry->index = _entries.size();
  _entries.push_back(std::move(entry));

  size_t sqIndex = _entries.back()->channel_group_id;
  assert(sqIndex >= _entries.front()->channel_group_id);
  sqIndex -= _entries.front()->channel_group_id;

  // As documented in deconvolutiontable.h, sqIndex may be at most 1 larger
  // than the largest existing sqIndex.
  assert(sqIndex <= _squaredGroups.size());
  if (sqIndex == _squaredGroups.size()) {
    // Create a new group for the entry.
    // The first entry of a group must have a psf image accessor.
    assert(_entries.back()->psf_accessor);
    _squaredGroups.emplace_back(1, _entries.back().get());
  } else {
    // Add the entry to an existing group.
    _squaredGroups[sqIndex].push_back(_entries.back().get());
  }
}
