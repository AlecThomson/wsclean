#include "griddingtaskfactory.h"

#include "../idg/averagebeam.h"
#include "../io/imagefilename.h"

GriddingTask GriddingTaskFactory::CreateBase(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    bool is_first_task) {
  GriddingTask task;

  task.observationInfo = observation_info_;
  task.isFirstTask = is_first_task;

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
  task.facetGroupIndex = entry.facetGroupIndex;

  task.msList = ms_helper_.InitializeMsList(entry);
  task.imageWeights = image_weight_initializer_.Initialize(entry, task.msList,
                                                           image_weight_cache);

  return task;
}

GriddingTask GriddingTaskFactory::CreatePsfTask(
    const ImagingTableEntry& entry, ImageWeightCache& image_weight_cache,
    bool is_first_task) {
  GriddingTask task = CreateBase(entry, image_weight_cache, is_first_task);
  task.operation = GriddingTask::Invert;
  task.imagePSF = true;
  task.polarization = entry.polarization;
  task.subtractModel = false;
  task.storeImagingWeights =
      image_weight_initializer_.GetSettings().writeImagingWeightSpectrumColumn;
  return task;
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
