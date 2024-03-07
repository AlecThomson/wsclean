#include "griddingtaskmanager.h"

#include <numeric>

#include "mpischeduler.h"
#include "threadedscheduler.h"

#include "../main/settings.h"

#include "../gridding/msgridderbase.h"
#include "../gridding/wsmsgridder.h"
#include "../gridding/directmsgridder.h"

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"

#include <schaapcommon/facets/facet.h>

#include "../wgridder/wgriddingmsgridder.h"

GriddingTaskManager::GriddingTaskManager(const Settings& settings)
    : _settings(settings), _writerLockManager(this) {}

GriddingTaskManager::~GriddingTaskManager() {}

std::unique_ptr<GriddingTaskManager> GriddingTaskManager::Make(
    const Settings& settings) {
  if (settings.UseMpi()) {
#ifdef HAVE_MPI
    return std::make_unique<MPIScheduler>(settings);
#else
    throw std::runtime_error("MPI not available");
#endif
  } else if (settings.parallelGridding > 1) {
    return std::make_unique<ThreadedScheduler>(settings);
  } else {
    return std::make_unique<GriddingTaskManager>(settings);
  }
}

Resources GriddingTaskManager::GetResources() const {
  return Resources(
      _settings.threadCount,
      GetAvailableMemory(_settings.memFraction, _settings.absMemLimit));
}

void GriddingTaskManager::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  std::vector<size_t> facet_indices(task.facets.size());
  std::iota(facet_indices.begin(), facet_indices.end(), 0);

  GriddingResult result;
  result.facets.resize(task.facets.size());
  std::mutex result_mutex;

  RunDirect(task, facet_indices, GetResources(), result, result_mutex);

  finishCallback(result);
}

void GriddingTaskManager::RunDirect(GriddingTask& task,
                                    const std::vector<size_t>& facet_indices,
                                    const Resources& resources,
                                    GriddingResult& result,
                                    std::mutex& result_mutex) {
  assert(!facet_indices.empty());
  assert(result.facets.size() == task.facets.size());

  std::unique_ptr<MSGridderBase> gridder;

  for (size_t facet_index : facet_indices) {
    assert(facet_index < task.facets.size());
    GriddingTask::FacetData& facet_task = task.facets[facet_index];
    GriddingResult::FacetData& facet_result = result.facets[facet_index];

    // Create a new gridder for each facet / sub-task, since gridders do not
    // support reusing them for multiple tasks.
    gridder = ConstructGridder(resources);
    InitializeGridderForTask(*gridder, task);

    const bool has_input_average_beam(facet_task.averageBeam);
    if (has_input_average_beam) {
      assert(dynamic_cast<IdgMsGridder*>(gridder.get()));
      IdgMsGridder& idgGridder = static_cast<IdgMsGridder&>(*gridder);
      idgGridder.SetAverageBeam(std::move(facet_task.averageBeam));
    }

    InitializeGridderForFacet(*gridder, facet_task);

    if (task.operation == GriddingTask::Invert) {
      gridder->Invert();
    } else {
      gridder->Predict(std::move(facet_task.modelImages));
    }

    // Add facet-specific result values to the result.
    facet_result.images = gridder->ResultImages();
    facet_result.actualWGridSize = gridder->ActualWGridSize();
    facet_result.averageCorrection = gridder->AverageCorrection();
    facet_result.averageH5Correction = gridder->AverageH5Correction();
    facet_result.cache = gridder->AcquireMetaDataCache();

    // The gridder resets visibility counters in each gridding invocation,
    // so they only contain the statistics of that invocation.
    facet_result.imageWeight = gridder->ImageWeight();
    facet_result.normalizationFactor = gridder->NormalizationFactor();
    facet_result.effectiveGriddedVisibilityCount =
        gridder->EffectiveGriddedVisibilityCount();
    {
      std::lock_guard<std::mutex> result_lock(result_mutex);
      result.griddedVisibilityCount += gridder->GriddedVisibilityCount();
      result.visibilityWeightSum += gridder->VisibilityWeightSum();
    }

    // If the average beam already exists on input, IDG will not recompute it,
    // so in that case there is no need to return the unchanged average beam.
    IdgMsGridder* idgGridder = dynamic_cast<IdgMsGridder*>(gridder.get());
    if (idgGridder && !has_input_average_beam) {
      facet_result.averageBeam = idgGridder->ReleaseAverageBeam();
    }
  }

  if (facet_indices.front() == 0) {
    // Store result values that are equal for all facets.
    result.unique_id = task.unique_id;
    result.startTime = gridder->StartTime();
    result.beamSize = gridder->BeamSize();
  }
}

std::unique_ptr<MSGridderBase> GriddingTaskManager::ConstructGridder(
    const Resources& resources) const {
  switch (_settings.gridderType) {
    case GridderType::IDG:
      return std::make_unique<IdgMsGridder>(_settings, resources);
    case GridderType::WGridder:
      return std::make_unique<WGriddingMSGridder>(_settings, resources, false);
    case GridderType::TunedWGridder:
      return std::make_unique<WGriddingMSGridder>(_settings, resources, true);
    case GridderType::DirectFT:
      switch (_settings.directFTPrecision) {
        case DirectFTPrecision::Float:
          return std::make_unique<DirectMSGridder<float>>(_settings, resources);
        case DirectFTPrecision::Double:
          return std::make_unique<DirectMSGridder<double>>(_settings,
                                                           resources);
        case DirectFTPrecision::LongDouble:
          return std::make_unique<DirectMSGridder<long double>>(_settings,
                                                                resources);
      }
      break;
    case GridderType::WStacking:
      return std::make_unique<WSMSGridder>(_settings, resources);
  }
  return {};
}

void GriddingTaskManager::InitializeGridderForTask(MSGridderBase& gridder,
                                                   const GriddingTask& task) {
  gridder.SetGridMode(_settings.gridMode);

  gridder.SetFacetGroupIndex(task.facetGroupIndex);
  gridder.SetImagePadding(_settings.imagePadding);
  gridder.SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder.SetPhaseCentreRA(task.observationInfo.phaseCentreRA);

  if (_settings.hasShift) {
    double main_image_dl = 0.0;
    double main_image_dm = 0.0;
    aocommon::ImageCoordinates::RaDecToLM(_settings.shiftRA, _settings.shiftDec,
                                          task.observationInfo.phaseCentreRA,
                                          task.observationInfo.phaseCentreDec,
                                          main_image_dl, main_image_dm);
    gridder.SetMainImageDL(main_image_dl);
    gridder.SetMainImageDM(main_image_dm);
  }

  gridder.SetPolarization(task.polarization);
  gridder.SetIsComplex(task.polarization == aocommon::Polarization::XY ||
                       task.polarization == aocommon::Polarization::YX);
  gridder.SetIsFirstTask(task.isFirstTask);
  gridder.SetImageWeights(task.imageWeights.get());
  if (task.operation == GriddingTask::Invert) {
    if (task.imagePSF) {
      if (_settings.ddPsfGridWidth > 1 || _settings.ddPsfGridHeight > 1) {
        gridder.SetPsfMode(PsfMode::kDirectionDependent);
      } else {
        gridder.SetPsfMode(PsfMode::kSingle);
      }
    } else {
      gridder.SetPsfMode(PsfMode::kNone);
    }
    gridder.SetDoSubtractModel(task.subtractModel);
    gridder.SetStoreImagingWeights(task.storeImagingWeights);
  } else {
    gridder.SetWriterLockManager(_writerLockManager);
  }

  gridder.ClearMeasurementSetList();
  for (const std::unique_ptr<MSDataDescription>& description : task.msList) {
    gridder.AddMeasurementSet(description->GetProvider(),
                              description->Selection());
  }
}

void GriddingTaskManager::InitializeGridderForFacet(
    MSGridderBase& gridder, GriddingTask::FacetData& facet_task) {
  const schaapcommon::facets::Facet* facet = facet_task.facet.get();
  gridder.SetIsFacet(facet != nullptr);
  if (facet) {
    gridder.SetFacetIndex(facet_task.index);
    gridder.SetImageWidth(facet->GetUntrimmedBoundingBox().Width());
    gridder.SetImageHeight(facet->GetUntrimmedBoundingBox().Height());
    gridder.SetTrimSize(facet->GetTrimmedBoundingBox().Width(),
                        facet->GetTrimmedBoundingBox().Height());
    gridder.SetFacetDirection(facet->RA(), facet->Dec());
  } else {
    gridder.SetImageWidth(_settings.paddedImageWidth);
    gridder.SetImageHeight(_settings.paddedImageHeight);
    gridder.SetTrimSize(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight);
  }
  gridder.SetLShift(facet_task.l_shift);
  gridder.SetMShift(facet_task.m_shift);

  std::unique_ptr<MetaDataCache> cache = std::move(facet_task.cache);
  if (!cache) cache = std::make_unique<MetaDataCache>();
  gridder.SetMetaDataCache(std::move(cache));
}