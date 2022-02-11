#include "deconvolutiontable.h"

#include <algorithm>
#include <cassert>

DeconvolutionTable::DeconvolutionTable(size_t n_original_groups,
                                       size_t n_deconvolution_groups,
                                       size_t channel_index_offset)
    : entries_(),
    channel_index_offset_(channel_index_offset),
      original_groups_(n_original_groups),
      deconvolution_groups_(
          (n_deconvolution_groups == 0)
              ? n_original_groups
              : std::min(n_original_groups, n_deconvolution_groups)) {}

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  const size_t original_channel_index = entry->original_channel_index;
  assert(original_channel_index < original_groups_.size());

  entry->index = entries_.size();
  entries_.push_back(std::move(entry));

  original_groups_[original_channel_index].push_back(entries_.back().get());

  const size_t deconvolution_index =
      (original_channel_index * deconvolution_groups_.size()) /
      original_groups_.size();
  deconvolution_groups_[deconvolution_index].push_back(entries_.back().get());
}
