#include "msgriddermanager.h"

#include <mutex>
#include <vector>

#include <aocommon/logger.h>
#include <schaapcommon/facets/facet.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include "directmsgridder.h"
#include "h5solutiondata.h"
#include "msgridderbase.h"
#include "msprovidercollection.h"
#include "wsmsgridder.h"

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"
#include "../main/settings.h"
#include "../structures/resources.h"
#include "../wgridder/wgriddingmsgridder.h"

using aocommon::Logger;

namespace wsclean {

void MSGridderManager::InitializeMS(GriddingTask& task) {
  for (const MsListItem& item : task.msList) {
    ms_provider_collection_.Add(item.ms_description->GetProvider(),
                                item.ms_description->Selection(),
                                item.ms_index);
  }
  ms_provider_collection_.InitializeMS();
}

void MSGridderManager::InitializeGridders(
    GriddingTask& task, const std::vector<size_t>& facet_indices,
    const Resources& resources,
    std::vector<GriddingResult::FacetData>& facet_results,
    GriddingTaskManager* writer_lock_manager) {
  for (size_t facet_index : facet_indices) {
    assert(facet_index < task.facets.size());

    // Create a new gridder for each facet / sub-task, since gridders do not
    // support reusing them for multiple tasks.
    std::unique_ptr<MSGridderBase> gridder = ConstructGridder(resources);
    GriddingTask::FacetData& facet_task = task.facets[facet_index];
    GriddingResult::FacetData& facet_result = facet_results[facet_index];

    if (solution_data_.HasData()) {
      gridder->SetH5Parm(
          solution_data_.GetH5Parms(), solution_data_.GetFirstSolutions(),
          solution_data_.GetSecondSolutions(), solution_data_.GetGainTypes());
    }
    InitializeGridderForTask(*gridder, task, writer_lock_manager);

    const bool has_input_average_beam(facet_task.averageBeam);
    if (has_input_average_beam) {
      assert(dynamic_cast<IdgMsGridder*>(gridder.get()));
      IdgMsGridder& idgGridder = static_cast<IdgMsGridder&>(*gridder);
      idgGridder.SetAverageBeam(std::move(facet_task.averageBeam));
    }

    InitializeGridderForFacet(*gridder, facet_task);

    facet_tasks_.emplace_back(
        GriddingFacetTask{std::move(gridder), facet_task, facet_result});
  }
}

void MSGridderManager::Invert() {
  InitializeMSDataVectors();

  for (const GriddingFacetTask& task : facet_tasks_) {
    task.facet_gridder->calculateOverallMetaData();
    task.facet_gridder->Invert();
  }
}

void MSGridderManager::Predict() {
  InitializeMSDataVectors();

  for (const GriddingFacetTask& task : facet_tasks_) {
    task.facet_gridder->calculateOverallMetaData();
    task.facet_gridder->Predict(std::move(task.facet_task.modelImages));
  }
}

void MSGridderManager::ProcessResults(std::mutex& result_mutex,
                                      GriddingResult& result,
                                      bool store_common_info) {
  for (auto& [gridder, facet_task, facet_result] : facet_tasks_) {
    // Add facet-specific result values to the reult.
    facet_result.images = gridder->ResultImages();
    facet_result.actualWGridSize = gridder->ActualWGridSize();
    facet_result.averageCorrection = gridder->GetAverageCorrection();
    facet_result.averageBeamCorrection = gridder->GetAverageBeamCorrection();
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
    const bool has_input_average_beam(facet_task.averageBeam);
    IdgMsGridder* idgGridder = dynamic_cast<IdgMsGridder*>(gridder.get());
    if (idgGridder && !has_input_average_beam) {
      facet_result.averageBeam = idgGridder->ReleaseAverageBeam();
    }

    if (store_common_info) {
      // Store result values that are equal for all facets.
      result.startTime = ms_provider_collection_.StartTime();
      result.beamSize = gridder->BeamSize();
    }
  }
}

std::unique_ptr<MSGridderBase> MSGridderManager::ConstructGridder(
    const Resources& resources) {
  switch (settings_.gridderType) {
    case GridderType::IDG:
      return std::make_unique<IdgMsGridder>(settings_, resources,
                                            ms_provider_collection_);
    case GridderType::WGridder:
      return std::make_unique<WGriddingMSGridder>(
          settings_, resources, ms_provider_collection_, false);
    case GridderType::TunedWGridder:
      return std::make_unique<WGriddingMSGridder>(
          settings_, resources, ms_provider_collection_, true);
    case GridderType::DirectFT:
      switch (settings_.directFTPrecision) {
        case DirectFTPrecision::Float:
          return std::make_unique<DirectMSGridder<float>>(
              settings_, resources, ms_provider_collection_);
        case DirectFTPrecision::Double:
          return std::make_unique<DirectMSGridder<double>>(
              settings_, resources, ms_provider_collection_);
        case DirectFTPrecision::LongDouble:
          return std::make_unique<DirectMSGridder<long double>>(
              settings_, resources, ms_provider_collection_);
      }
      break;
    case GridderType::WStacking:
      return std::make_unique<WSMSGridder>(settings_, resources,
                                           ms_provider_collection_);
  }
  return {};
}

void MSGridderManager::InitializeGridderForTask(
    MSGridderBase& gridder, const GriddingTask& task,
    GriddingTaskManager* writer_lock_manager) {
  gridder.SetGridMode(settings_.gridMode);

  gridder.SetFacetGroupIndex(task.facetGroupIndex);
  gridder.SetImagePadding(settings_.imagePadding);
  gridder.SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder.SetPhaseCentreRA(task.observationInfo.phaseCentreRA);

  if (settings_.hasShift) {
    double main_image_dl = 0.0;
    double main_image_dm = 0.0;
    aocommon::ImageCoordinates::RaDecToLM(settings_.shiftRA, settings_.shiftDec,
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
      if (settings_.ddPsfGridWidth > 1 || settings_.ddPsfGridHeight > 1) {
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
    gridder.SetWriterLockManager(writer_lock_manager);
  }
}

void MSGridderManager::InitializeGridderForFacet(
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
    gridder.SetImageWidth(settings_.paddedImageWidth);
    gridder.SetImageHeight(settings_.paddedImageHeight);
    gridder.SetTrimSize(settings_.trimmedImageWidth,
                        settings_.trimmedImageHeight);
  }
  gridder.SetLShift(facet_task.l_shift);
  gridder.SetMShift(facet_task.m_shift);

  std::unique_ptr<MetaDataCache> cache = std::move(facet_task.cache);
  if (!cache) cache = std::make_unique<MetaDataCache>();
  gridder.SetMetaDataCache(std::move(cache));
}

}  // namespace wsclean
