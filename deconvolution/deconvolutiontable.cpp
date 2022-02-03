#include "deconvolutiontable.h"

#include <cassert>

void DeconvolutionTable::AddEntry(
    std::unique_ptr<DeconvolutionTableEntry> entry) {
  entry->index = _entries.size();
  _entries.push_back(std::move(entry));

  size_t sqIndex = _entries.back()->squaredDeconvolutionIndex;
  assert(sqIndex >= _entries.front()->squaredDeconvolutionIndex);
  sqIndex -= _entries.front()->squaredDeconvolutionIndex;
  assert(sqIndex <= _squaredGroups.size());
  if (sqIndex == _squaredGroups.size()) {
    _squaredGroups.emplace_back(1, _entries.back());
  } else {
    _squaredGroups[sqIndex].push_back(_entries.back());
  }
}
