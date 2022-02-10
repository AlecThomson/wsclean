#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  const size_t channel_index = entry->channel_index;
  assert(channel_index < channel_groups_.size());

  entry->index = entries_.size();
  entries_.push_back(std::move(entry));

  channel_groups_[channel_index].push_back(entries_.back().get());
}
