#include "mshelper.h"

#include <mutex>

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>
#include <aocommon/logger.h>

using aocommon::Logger;

void MsHelper::PerformReordering(const ImagingTable& imaging_table,
                                 bool is_predict_mode) {
  std::mutex mutex;

  assert(reordered_ms_handles_.empty());

  reordered_ms_handles_.resize(settings_.filenames.size());
  bool use_model = settings_.deconvolutionMGain != 1.0 || is_predict_mode ||
                   settings_.subtractModel || settings_.continuedRun;
  bool initial_model_required =
      settings_.subtractModel || settings_.continuedRun;

  if (settings_.parallelReordering != 1) Logger::Info << "Reordering...\n";

  aocommon::CountingSemaphore semaphore(settings_.parallelReordering);
  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, settings_.filenames.size(), [&](size_t ms_index) {
    const aocommon::MultiBandData& band_data = ms_bands_[ms_index];
    aocommon::ScopedCountingSemaphoreLock semaphore_lock(semaphore);
    std::vector<ReorderedMs::ChannelRange> channels;
    // The partIndex needs to increase per data desc ids and channel ranges
    std::map<aocommon::PolarizationEnum, size_t> next_index;
    for (size_t sq_index = 0; sq_index != imaging_table.SquaredGroupCount();
         ++sq_index) {
      const ImagingTable::Groups facet_groups =
          imaging_table.FacetGroups([sq_index](const ImagingTableEntry& e) {
            return e.squaredDeconvolutionIndex == sq_index;
          });
      for (const ImagingTable::Group& facet_group : facet_groups) {
        // The band information is determined from the first facet in the group.
        // After this, all facet entries inside the group are updated.
        const ImagingTableEntry& entry = *facet_group.front();
        for (size_t d = 0; d != band_data.DataDescCount(); ++d) {
          MSSelection selection{global_selection_};
          const size_t band_index = band_data.GetBandIndex(d);

          if (settings_.IsBandSelected(band_index) &&
              selection.SelectMsChannels(band_data, d, entry)) {
            if (entry.polarization == *settings_.polarizations.begin()) {
              ReorderedMs::ChannelRange r;
              r.dataDescId = d;
              r.start = selection.ChannelRangeStart();
              r.end = selection.ChannelRangeEnd();
              channels.push_back(r);
            }
            for (const std::shared_ptr<ImagingTableEntry>& facet_entry :
                 facet_group) {
              facet_entry->msData[ms_index].bands[d].partIndex =
                  next_index[entry.polarization];
            }
            ++next_index[entry.polarization];
          }
        }
      }
    }

    ReorderedMs::Handle part_ms = ReorderedMs::Partition(
        settings_.filenames[ms_index], channels, global_selection_,
        settings_.dataColumnName, use_model, initial_model_required, settings_);
    std::lock_guard<std::mutex> lock(mutex);
    reordered_ms_handles_[ms_index] = std::move(part_ms);
    if (settings_.parallelReordering != 1)
      Logger::Info << "Finished reordering " << settings_.filenames[ms_index]
                   << " [" << ms_index << "]\n";
  });
}

std::vector<std::unique_ptr<MSDataDescription>> MsHelper::InitializeMsList(
    const ImagingTableEntry& entry) const {
  const aocommon::PolarizationEnum polarization =
      settings_.GetProviderPolarization(entry.polarization);

  std::vector<std::unique_ptr<MSDataDescription>> ms_list;

  for (size_t ms_index = 0; ms_index != settings_.filenames.size();
       ++ms_index) {
    const aocommon::MultiBandData& band_data = ms_bands_[ms_index];

    for (size_t data_description_id = 0;
         data_description_id != band_data.DataDescCount();
         ++data_description_id) {
      MSSelection selection{global_selection_};
      const size_t band_index = band_data.GetBandIndex(data_description_id);

      if (settings_.IsBandSelected(band_index) &&
          selection.SelectMsChannels(band_data, data_description_id, entry)) {
        std::unique_ptr<MSDataDescription> data_description;
        if (settings_.doReorder)
          data_description = MSDataDescription::ForReordered(
              reordered_ms_handles_[ms_index], selection,
              entry.msData[ms_index].bands[data_description_id].partIndex,
              polarization, data_description_id, settings_.UseMpi());
        else
          data_description = MSDataDescription::ForContiguous(
              settings_.filenames[ms_index], settings_.dataColumnName,
              selection, polarization, data_description_id, settings_.UseMpi());
        ms_list.emplace_back(std::move(data_description));
      }
    }
  }

  return ms_list;
}
