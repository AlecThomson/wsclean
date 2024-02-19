#include "griddingtaskmanager.h"

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
  GriddingResult result = RunDirect(std::move(task), GetResources());
  finishCallback(result);
}

GriddingResult GriddingTaskManager::RunDirect(GriddingTask&& task,
                                              const Resources& resources) {
  std::unique_ptr<MSGridderBase> gridder = ConstructGridder(resources);
  gridder->SetGridMode(_settings.gridMode);

  GriddingTask::FacetData& facet_task = task.facets.front();

  gridder->ClearMeasurementSetList();
  std::vector<std::unique_ptr<MSProvider>> msProviders;
  for (auto& p : task.msList) {
    msProviders.emplace_back(p->GetProvider());
    gridder->AddMeasurementSet(msProviders.back().get(), p->Selection());
  }

  const bool has_input_average_beam(facet_task.averageBeam);
  if (has_input_average_beam) {
    assert(dynamic_cast<IdgMsGridder*>(gridder.get()));
    IdgMsGridder& idgGridder = static_cast<IdgMsGridder&>(*gridder);
    idgGridder.SetAverageBeam(std::move(facet_task.averageBeam));
  }

  gridder->SetFacetGroupIndex(task.facetGroupIndex);
  const schaapcommon::facets::Facet* facet = facet_task.facet.get();
  gridder->SetIsFacet(facet != nullptr);
  if (facet) {
    gridder->SetFacetIndex(facet_task.index);
    gridder->SetImageWidth(facet->GetUntrimmedBoundingBox().Width());
    gridder->SetImageHeight(facet->GetUntrimmedBoundingBox().Height());
    gridder->SetTrimSize(facet->GetTrimmedBoundingBox().Width(),
                         facet->GetTrimmedBoundingBox().Height());
    gridder->SetFacetDirection(facet->RA(), facet->Dec());
  } else {
    gridder->SetImageWidth(_settings.paddedImageWidth);
    gridder->SetImageHeight(_settings.paddedImageHeight);
    gridder->SetTrimSize(_settings.trimmedImageWidth,
                         _settings.trimmedImageHeight);
  }
  gridder->SetLShift(facet_task.l_shift);
  gridder->SetMShift(facet_task.m_shift);
  std::unique_ptr<MetaDataCache> cache = std::move(facet_task.cache);
  if (!cache) cache = std::make_unique<MetaDataCache>();
  gridder->SetMetaDataCache(std::move(cache));

  gridder->SetImagePadding(_settings.imagePadding);
  gridder->SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder->SetPhaseCentreRA(task.observationInfo.phaseCentreRA);

  if (_settings.hasShift) {
    double main_image_dl = 0.0;
    double main_image_dm = 0.0;
    aocommon::ImageCoordinates::RaDecToLM(_settings.shiftRA, _settings.shiftDec,
                                          task.observationInfo.phaseCentreRA,
                                          task.observationInfo.phaseCentreDec,
                                          main_image_dl, main_image_dm);
    gridder->SetMainImageDL(main_image_dl);
    gridder->SetMainImageDM(main_image_dm);
  }

  gridder->SetPolarization(task.polarization);
  gridder->SetIsComplex(task.polarization == aocommon::Polarization::XY ||
                        task.polarization == aocommon::Polarization::YX);
  gridder->SetIsFirstTask(task.isFirstTask);
  gridder->SetImageWeights(task.imageWeights.get());
  if (task.operation == GriddingTask::Invert) {
    if (task.imagePSF) {
      if (_settings.ddPsfGridWidth > 1 || _settings.ddPsfGridHeight > 1) {
        gridder->SetPsfMode(PsfMode::kDirectionDependent);
      } else {
        gridder->SetPsfMode(PsfMode::kSingle);
      }
    } else {
      gridder->SetPsfMode(PsfMode::kNone);
    }
    gridder->SetDoSubtractModel(task.subtractModel);
    gridder->SetStoreImagingWeights(task.storeImagingWeights);
    gridder->Invert();
  } else {
    gridder->SetWriterLockManager(_writerLockManager);
    gridder->Predict(std::move(facet_task.modelImages));
  }

  GriddingResult result;
  GriddingResult::FacetData& facet_result = result.facets.front();
  facet_result.images = gridder->ResultImages();
  result.unique_id = task.unique_id;
  result.startTime = gridder->StartTime();
  result.beamSize = gridder->BeamSize();
  facet_result.imageWeight = gridder->ImageWeight();
  facet_result.normalizationFactor = gridder->NormalizationFactor();
  facet_result.actualWGridSize = gridder->ActualWGridSize();
  result.griddedVisibilityCount = gridder->GriddedVisibilityCount();
  facet_result.effectiveGriddedVisibilityCount =
      gridder->EffectiveGriddedVisibilityCount();
  result.visibilityWeightSum = gridder->VisibilityWeightSum();
  facet_result.averageCorrection = gridder->AverageCorrection();
  facet_result.averageH5Correction = gridder->AverageH5Correction();
  facet_result.cache = gridder->AcquireMetaDataCache();

  // If the average beam already exists on input, IDG will not recompute it, so
  // in that case there is no need to return the unchanged average beam.
  IdgMsGridder* idgGridder = dynamic_cast<IdgMsGridder*>(gridder.get());
  if (idgGridder && !has_input_average_beam) {
    facet_result.averageBeam = idgGridder->ReleaseAverageBeam();
  }
  return result;
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
