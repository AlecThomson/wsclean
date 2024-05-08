
#include "wgriddingmsgridder.h"

#include "wgriddinggridder_simple.h"

#include "../msproviders/msreaders/msreader.h"

#include "../msproviders/msprovider.h"

#include "../system/buffered_lane.h"

#include "../structures/imageweights.h"

#include "BS_thread_pool.hpp"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <schaapcommon/fft/resampler.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <algorithm>

using aocommon::Image;
using aocommon::Logger;

WGriddingMSGridder::WGriddingMSGridder(const Settings& settings,
                                       const Resources& resources,
                                       bool use_tuned_wgridder)
    : MSGridderBase(settings),
      resources_(resources),
      accuracy_(GetSettings().wgridderAccuracy),
      use_tuned_wgridder_(use_tuned_wgridder) {
  // It may happen that several schaapcommon::fft::Resamplers are created
  // concurrently, so we must make sure that the FFTW planner can deal with
  // this.
  fftwf_make_planner_thread_safe();
}

WGriddingMSGridder::~WGriddingMSGridder() = default;

std::unique_ptr<WGriddingGridderBase> WGriddingMSGridder::MakeGridder(
    size_t width, size_t height) const {
  if (accuracy_ <= 1.01e-5) {
    return std::make_unique<WGriddingGridder_Simple<double>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        resources_.NCpus() / num_parallel_nested_, accuracy_, 0,
        use_tuned_wgridder_);
  } else {
    return std::make_unique<WGriddingGridder_Simple<float>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        resources_.NCpus() / num_parallel_nested_, accuracy_, 0,
        use_tuned_wgridder_);
  }
}

size_t WGriddingMSGridder::calculateMaxNRowsInMemory(
    size_t channelCount) const {
  size_t constantMem;
  size_t perVisMem;
  gridder_->MemoryUsage(constantMem, perVisMem);
  if (int64_t(constantMem) >= resources_.Memory()) {
    // Assume that half the memory is necessary for the constant parts (like
    // image grid), and the other half remains available for the dynamic buffers
    constantMem = resources_.Memory() / 2;
    Logger::Warn << "Not enough memory available for doing the gridding:\n"
                    "swapping might occur!\n";
  }
  const uint64_t memForBuffers = resources_.Memory() - constantMem;

  const uint64_t memPerRow =
      perVisMem                                       // ducc internal usage
      + (sizeof(std::complex<float>) * channelCount)  // visibilities
      + (sizeof(float) * channelCount)                // weights
      + sizeof(double) * 3                            // uvw
      + sizeof(MSProvider::MetaData);                 // metadata
  const size_t maxNRows = std::max(memForBuffers / memPerRow, uint64_t(100));
  if (maxNRows < 1000) {
    Logger::Warn << "Less than 1000 data rows fit in memory: this probably "
                    "means performance is going to be very poor!\n";
  }

  return maxNRows;
}

void WGriddingMSGridder::gridMeasurementSet(MSData& msData) {
  const size_t num_gridding_threads = std::min(
      GetSettings().parallelGridding, static_cast<size_t>(numNestedGridders()));
  const size_t num_threads_per_gridder =
      GetSettings().threadCount / GetSettings().parallelGridding;
  const size_t num_overall_threads =
      num_gridding_threads * num_threads_per_gridder;

  BS::thread_pool all_core_thread_pool(num_overall_threads);
  BS::thread_pool gridder_task_thread_grid_pool(num_gridding_threads);

  Logger::Info << "gridMeasurementSet..\n";
  const aocommon::BandData selectedBand(msData.SelectedBand());
  for (auto& gridder : nestedGridders) {
    gridder->StartMeasurementSet(msData, false);
  }

  const size_t n_vis_polarizations = msData.ms_provider->NPolarizations();
  const size_t data_size = selectedBand.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> modelBuffer(data_size);
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  size_t totalNRows = 0;
  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  aocommon::UVector<double> uvw_buffer(maxNRows * 3);
  aocommon::UVector<std::complex<float>> visibility_buffer(maxNRows *
                                                           data_size);
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<MSProvider::MetaData> metadata_buffer(maxNRows);

  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();

  // Iterate over chunks until all data has been gridded
  while (msReader->CurrentRowAvailable()) {
    Logger::Debug << "Max " << maxNRows << " rows fit in memory.\n";
    Logger::Info << "Loading data in memory...\n";

    size_t nRows = 0;
    std::complex<float>* vis_buffer_index = visibility_buffer.data();
    double* uvw_buffer_index = uvw_buffer.data();
    MSProvider::MetaData* meta_data = metadata_buffer.data();

    // Read / fill the chunk
    while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
      msReader->ReadMeta(*meta_data);
      uvw_buffer_index[0] = meta_data->uInM;
      uvw_buffer_index[1] = meta_data->vInM;
      uvw_buffer_index[2] = meta_data->wInM;

      // Read and store all visibilities and weights, we need them all in memory
      // when calling 'InlineApplyWeightsAndCorrections' so that we can
      // calculate the final visibilities to return correctly
      ReadVisibilities(*msReader, vis_buffer_index, weight_buffer.data(),
                       modelBuffer.data());
      CalculateWeights(uvw_buffer_index, vis_buffer_index, selectedBand,
                       weight_buffer.data(), modelBuffer.data(),
                       isSelected.data(), *meta_data, msData.antenna_names,
                       nestedGridders);

      // We don't currently store the imaging weights
      // For testing we don't need this
      // Final implementation will need this
      // if (StoreImagingWeights())
      // ms_reader.WriteImagingWeights(scratch_image_weights_.data());

      vis_buffer_index += data_size;
      uvw_buffer_index += 3;
      ++meta_data;

      ++nRows;
      msReader->NextInputRow();
    }

    Logger::Info << "Gridding " << nRows << " rows for "
                 << nestedGridders.size() << " facets across "
                 << num_gridding_threads << "threads...\n";

    for (size_t facet_index = 0; facet_index < nestedGridders.size();
         ++facet_index) {
      gridder_task_thread_grid_pool.detach_task([&, facet_index] {
        std::stringstream logMessage;
        logMessage << "Gridding facet" << facet_index << "\n";
        Logger::Info << logMessage.str();
        const auto& gridder = nestedGridders[facet_index];

        gridder->gridded_visibility_count_ = gridded_visibility_count_;
        gridder->visibility_weight_sum_ = visibility_weight_sum_;
        gridder->max_gridded_weight_ = max_gridded_weight_;
        gridder->total_weight_ = total_weight_;

        // Create a callback and virtual memory buffer object that can apply
        // weights and corrections on the fly when passed into DUCC
        std::function<std::complex<float>&(size_t index)>
            getVisibilityWithWeight;
        if (n_vis_polarizations == 1) {
          getVisibilityWithWeight = [&](size_t index) -> std::complex<float>& {
            size_t row = index / data_size;
            size_t channel = index % data_size;
            MSProvider::MetaData& meta_data = metadata_buffer[row];

            return gridder->InlineApplyWeightsAndCorrections<1>(
                msData.antenna_names, selectedBand, visibility_buffer.data(),
                index, row, channel, meta_data);
          };
        }
        auto ms =
            cmav<std::complex<float>, 2, cvirtmembuf<std::complex<float>>>(
                {nRows, selectedBand.ChannelCount()},
                cvirtmembuf<std::complex<float>>(getVisibilityWithWeight));

        ((WGriddingMSGridder*)gridder)
            ->gridder_->AddInversionData(nRows, selectedBand.ChannelCount(),
                                         uvw_buffer.data(), frequencies.data(),
                                         ms);
        Logger::Info << "Done gridding facet" << facet_index << "\n";
      });
    }
    gridder_task_thread_grid_pool.wait();

    // Clear 'per chunk' data that we have cached so that they are fresh for the
    // next chunk of data
    for (size_t facet_index = 0; facet_index < nestedGridders.size();
         ++facet_index) {
      const auto& gridder = nestedGridders[facet_index];
      gridder->time_offsets_.clear();
    }

    totalNRows += nRows;
  }  // end of chunk

  msData.totalRowsProcessed += totalNRows;
}

void WGriddingMSGridder::predictMeasurementSet(MSData& msData) {
  msData.ms_provider->ReopenRW();
  const aocommon::BandData selectedBand(msData.SelectedBand());
  StartMeasurementSet(msData, true);

  size_t totalNRows = 0;

  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  aocommon::UVector<double> uvwBuffer(maxNRows * 3);
  // Iterate over chunks until all data has been gridded
  msData.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();
  while (msReader->CurrentRowAvailable()) {
    size_t nRows = 0;

    // Read / fill the chunk
    Logger::Info << "Loading metadata...\n";
    // Read from metadata buffer
    std::vector<MSProvider::MetaData> metaDataBuffer;
    while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
      metaDataBuffer.emplace_back();
      ReadPredictMetaData(metaDataBuffer[nRows]);
      uvwBuffer[nRows * 3] = metaDataBuffer[nRows].uInM;
      uvwBuffer[nRows * 3 + 1] = metaDataBuffer[nRows].vInM;
      uvwBuffer[nRows * 3 + 2] = metaDataBuffer[nRows].wInM;
      nRows++;

      msReader->NextInputRow();
    }

    Logger::Info << "Predicting " << nRows << " rows...\n";
    aocommon::UVector<std::complex<float>> visBuffer(
        maxNRows * selectedBand.ChannelCount());
    gridder_->PredictVisibilities(nRows, selectedBand.ChannelCount(),
                                  uvwBuffer.data(), frequencies.data(),
                                  visBuffer.data());

    Logger::Info << "Writing...\n";
    for (size_t row = 0; row != nRows; ++row) {
      WriteCollapsedVisibilities(
          *msData.ms_provider, msData.antenna_names, selectedBand,
          &visBuffer[row * selectedBand.ChannelCount()], metaDataBuffer[row]);
    }
    totalNRows += nRows;
  }  // end of chunk

  msData.totalRowsProcessed += totalNRows;
}

void WGriddingMSGridder::getActualTrimmedSize(size_t& trimmedWidth,
                                              size_t& trimmedHeight) const {
  trimmedWidth = std::ceil(ActualInversionWidth() / ImagePadding());
  trimmedHeight = std::ceil(ActualInversionHeight() / ImagePadding());

  // In facet-based imaging, the alignment is 4, see wsclean.cpp. Also for
  // monolithic imaging - in which just an even number would suffice -
  // the trimmedWidth and trimmedHeight are defined to be divisable by 4.
  const size_t alignment = 4;
  if (trimmedWidth % alignment != 0) {
    trimmedWidth += alignment - (trimmedWidth % alignment);
  }
  if (trimmedHeight % alignment != 0) {
    trimmedHeight += alignment - (trimmedHeight % alignment);
  }
  trimmedWidth = std::min(trimmedWidth, ActualInversionWidth());
  trimmedHeight = std::min(trimmedHeight, ActualInversionHeight());
}

void WGriddingMSGridder::Invert() {
  std::vector<std::vector<MSData>> msDataVectors;
  for (const auto& pgridder : nestedGridders) {
    std::vector<MSData> msDataVector;
    pgridder->initializeMSDataVector(msDataVector);
    msDataVectors.emplace_back(msDataVector);
  }

  const size_t num_overall_threads = resources_.NCpus();
  const size_t num_gridding_threads = num_parallel_nested_;

  BS::thread_pool all_core_thread_pool(num_overall_threads);
  for (size_t gridder_index = 0; gridder_index < nestedGridders.size();
       ++gridder_index) {
    all_core_thread_pool.detach_task([&, gridder_index] {
      WGriddingMSGridder* pgridder =
          (WGriddingMSGridder*)(nestedGridders[gridder_index]);
      auto msDataVector = msDataVectors[gridder_index];

      pgridder->calculateOverallMetaData(msDataVector);

      size_t trimmedWidth, trimmedHeight;
      pgridder->getActualTrimmedSize(trimmedWidth, trimmedHeight);

      pgridder->gridder_ = pgridder->MakeGridder(trimmedWidth, trimmedHeight);
      pgridder->gridder_->InitializeInversion();

      pgridder->resetVisibilityCounters();
    });
  }
  all_core_thread_pool.wait();

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    MSData& msData = msDataVectors[0][i];
    gridMeasurementSet(msData);
  }

  for (size_t gridder_index = 0; gridder_index < nestedGridders.size();
       ++gridder_index) {
    all_core_thread_pool.detach_task([&, gridder_index] {
      WGriddingMSGridder* pgridder =
          (WGriddingMSGridder*)(nestedGridders[gridder_index]);

      pgridder->gridder_->FinalizeImage(1.0 / pgridder->totalWeight());

      Logger::Info << "Gridded visibility count: "
                   << double(pgridder->GriddedVisibilityCount());
      if (pgridder->Weighting().IsNatural())
        Logger::Info << ", effective count after weighting: "
                     << pgridder->EffectiveGriddedVisibilityCount();
      Logger::Info << '\n';

      pgridder->_image = Image(pgridder->ActualInversionWidth(),
                               pgridder->ActualInversionHeight());
      {
        std::vector<float> imageFloat = pgridder->gridder_->RealImage();
        for (size_t i = 0; i < imageFloat.size(); ++i)
          pgridder->_image[i] = imageFloat[i];
      }

      if (pgridder->ImageWidth() != pgridder->ActualInversionWidth() ||
          pgridder->ImageHeight() != pgridder->ActualInversionHeight()) {
        // Interpolate the image
        // The input is of size ActualInversionWidth() x ActualInversionHeight()
        schaapcommon::fft::Resampler resampler(
            pgridder->ActualInversionWidth(), pgridder->ActualInversionHeight(),
            pgridder->ImageWidth(), pgridder->ImageHeight(),
            resources_.NCpus() / numNestedGridders());

        Image resized(pgridder->ImageWidth(), pgridder->ImageHeight());
        resampler.Resample(pgridder->_image.Data(), resized.Data());
        pgridder->_image = std::move(resized);
      }

      if (pgridder->TrimWidth() != pgridder->ImageWidth() ||
          pgridder->TrimHeight() != pgridder->ImageHeight()) {
        Logger::Debug << "Trimming " << pgridder->ImageWidth() << " x "
                      << pgridder->ImageHeight() << " -> "
                      << pgridder->TrimWidth() << " x "
                      << pgridder->TrimHeight() << '\n';

        pgridder->_image = pgridder->_image.Trim(pgridder->TrimWidth(),
                                                 pgridder->TrimHeight());
      }
    });
  }
  all_core_thread_pool.wait();
}

void WGriddingMSGridder::Predict(std::vector<Image>&& images) {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);
  calculateOverallMetaData(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getActualTrimmedSize(trimmedWidth, trimmedHeight);

  gridder_ = MakeGridder(trimmedWidth, trimmedHeight);

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmedImage(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    Image::Untrim(untrimmedImage.Data(), ImageWidth(), ImageHeight(),
                  images[0].Data(), TrimWidth(), TrimHeight());
    images[0] = std::move(untrimmedImage);
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    Image resampledImage(ImageWidth(), ImageHeight());
    schaapcommon::fft::Resampler resampler(
        ImageWidth(), ImageHeight(), ActualInversionWidth(),
        ActualInversionHeight(), resources_.NCpus());

    resampler.Resample(images[0].Data(), resampledImage.Data());
    images[0] = std::move(resampledImage);
  }

  gridder_->InitializePrediction(images[0].Data());
  images[0].Reset();

  for (MSData& msData : msDataVector) {
    predictMeasurementSet(msData);
  }
}
