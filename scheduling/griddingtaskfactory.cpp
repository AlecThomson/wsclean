#include "griddingtaskfactory.h"

#include "../idg/averagebeam.h"
#include "../io/imagefilename.h"

void GriddingTaskFactory::AddFacet(GriddingTask& task,
                                   const ImagingTableEntry& entry) {
  double l_shift = l_shift_;
  double m_shift = m_shift_;

  if (entry.facet) {
    const Settings& settings = image_weight_initializer_.GetSettings();
    l_shift -= entry.centreShiftX * settings.pixelScaleX;
    m_shift += entry.centreShiftY * settings.pixelScaleY;
  } else {
    assert(entry.facetIndex == 0);
  }
  task.facets.emplace_back(entry.facetIndex, l_shift, m_shift,
                           std::move(meta_data_cache_[entry.index]),
                           entry.facet);
}

GriddingTask GriddingTaskFactory::CreateBase(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    bool is_first_task) {
  GriddingTask task;

  task.observationInfo = observation_info_;
  task.isFirstTask = is_first_task;

  AddFacet(task, entry);
  task.facetGroupIndex = entry.facetGroupIndex;

  task.msList = ms_helper_.InitializeMsList(entry);
  task.imageWeights = image_weight_initializer_.Initialize(entry, task.msList,
                                                           image_weight_cache);

  return task;
}

std::vector<GriddingTask> GriddingTaskFactory::CreatePsfTasks(
    const ImagingTable::Group& facet_group,
    ImageWeightCache& image_weight_cache, bool combine_facets,
    bool is_first_task) {
  const bool store_imaging_weights =
      image_weight_initializer_.GetSettings().writeImagingWeightSpectrumColumn;

  std::vector<GriddingTask> tasks;
  tasks.reserve(combine_facets ? 1 : facet_group.size());

  for (const std::shared_ptr<ImagingTableEntry>& entry : facet_group) {
    // During PSF imaging, the average beam will never exist, so it is not
    // necessary to set the average beam in the task.

    if (combine_facets && !tasks.empty()) {
      assert(!tasks.front().facets.empty());
      assert(tasks.front().facets.front().facet);
      assert(entry->facet);
      assert(tasks.front().facetGroupIndex == entry->facetGroupIndex);
      AddFacet(tasks.front(), *entry);
    } else {  // Create a new task.
      tasks.push_back(CreateBase(*entry, image_weight_cache, is_first_task));
      tasks.back().operation = GriddingTask::Invert;
      tasks.back().imagePSF = true;
      tasks.back().polarization = entry->polarization;
      tasks.back().subtractModel = false;
      tasks.back().storeImagingWeights = store_imaging_weights;
    }

    is_first_task = false;
  }

  return tasks;
}

GriddingTask GriddingTaskFactory::CreateInvertTask(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    bool is_first_task, bool is_first_inversion,
    std::unique_ptr<AverageBeam> average_beam) {
  const Settings& settings = image_weight_initializer_.GetSettings();

  GriddingTask task = CreateBase(entry, image_weight_cache, is_first_task);

  task.operation = GriddingTask::Invert;
  task.imagePSF = false;
  if (settings.gridderType == GridderType::IDG &&
      settings.polarizations.size() != 1)
    task.polarization = aocommon::Polarization::FullStokes;
  else
    task.polarization = entry.polarization;
  task.subtractModel =
      !is_first_inversion || settings.subtractModel || settings.continuedRun;
  task.storeImagingWeights =
      is_first_inversion && settings.writeImagingWeightSpectrumColumn;

  GriddingTask::FacetData& facet_task = task.facets.front();
  facet_task.averageBeam = std::move(average_beam);

  return task;
}

GriddingTask GriddingTaskFactory::CreatePredictTask(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    std::vector<aocommon::Image>&& model_images,
    std::unique_ptr<AverageBeam> average_beam) {
  const Settings& settings = image_weight_initializer_.GetSettings();

  GriddingTask task = CreateBase(entry, image_weight_cache, false);
  task.operation = GriddingTask::Predict;

  if (settings.gridderType == GridderType::IDG &&
      settings.polarizations.size() != 1)
    task.polarization = aocommon::Polarization::FullStokes;
  else
    task.polarization = entry.polarization;

  task.storeImagingWeights = false;

  GriddingTask::FacetData& facet_task = task.facets.front();
  facet_task.modelImages = std::move(model_images);
  facet_task.averageBeam = std::move(average_beam);

  return task;
}
