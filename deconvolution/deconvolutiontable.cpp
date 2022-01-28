#include "deconvolutiontable.h"

#include <map>

void DeconvolutionTable::Update() {
  std::map<size_t, Group> group_map;

  for (const EntryPtr& e : _entries) {
    group_map[e->squaredDeconvolutionIndex].push_back(e);
  }

  _squaredGroups.clear();
  _squaredGroups.reserve(group_map.size());
  for (auto& item : group_map) {
    _squaredGroups.emplace_back(std::move(item.second));
  }
}
