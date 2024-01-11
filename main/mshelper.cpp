#include "mshelper.h"

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
          data_description = MSDataDescription::ForPartitioned(
              partitioned_ms_handles_[ms_index], selection,
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
