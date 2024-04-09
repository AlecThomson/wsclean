
#include "wgriddingmsgridder.h"

#include "wgriddinggridder_simple.h"

#include "../msproviders/msreaders/msreader.h"

#include "../msproviders/msprovider.h"

#include "../system/buffered_lane.h"

#include "../structures/imageweights.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <schaapcommon/fft/resampler.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

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
        resources_.NCpus(), accuracy_, 0, use_tuned_wgridder_);
  } else {
    return std::make_unique<WGriddingGridder_Simple<float>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        resources_.NCpus(), accuracy_, 0, use_tuned_wgridder_);
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

  const uint64_t memPerRow = (perVisMem + sizeof(std::complex<float>)) *
                                 channelCount       // vis themselves
                             + sizeof(double) * 3;  // uvw
  const size_t maxNRows = std::max(memForBuffers / memPerRow, uint64_t(100));
  if (maxNRows < 1000) {
    Logger::Warn << "Less than 1000 data rows fit in memory: this probably "
                    "means performance is going to be very poor!\n";
  }

  return maxNRows;
}

void WGriddingMSGridder::gridMeasurementSet(MSData& msData) {
  const aocommon::BandData selectedBand(msData.SelectedBand());
  StartMeasurementSet(msData, false);

  const size_t n_vis_polarizations = msData.ms_provider->NPolarizations();
  const size_t data_size = selectedBand.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> modelBuffer(data_size);
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  size_t totalNRows = 0;
  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  // For the sake of simplicity the prototype does not *yet* support more than 1
  // polarization here. We will initially be testing the prototype in a way that
  // n_vis_polarizations is always 1 We force an exit here to  make it apparent
  // if a mistake has been made As this develops into a full feature we will
  // have to change and implement this
  if (n_vis_polarizations > 1) {
    exit(657);
  }

  aocommon::UVector<double> uvw_buffer(maxNRows * 3);
  aocommon::UVector<std::complex<float>> visibility_buffer(maxNRows *
                                                           data_size);
  aocommon::UVector<float> weight_buffer(maxNRows * data_size);
  aocommon::UVector<float> uv_buffer(maxNRows * selectedBand.ChannelCount());
  aocommon::UVector<MSProvider::MetaData> metadata_buffer(maxNRows);

  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();

  // Iterate over chunks until all data has been gridded
  while (msReader->CurrentRowAvailable()) {
    Logger::Debug << "Max " << maxNRows << " rows fit in memory.\n";
    Logger::Info << "Loading data in memory...\n";

    size_t nRows = 0;
    std::complex<float>* vis_buffer_index = visibility_buffer.data();
    float* weight_buffer_index = weight_buffer.data();
    float* uv_buffer_index = uv_buffer.data();
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
      ReadVisibilities(*msReader, vis_buffer_index, weight_buffer_index,
                       modelBuffer.data());
      CalculateWeights(uvw_buffer_index, vis_buffer_index, selectedBand,
                       weight_buffer_index, modelBuffer.data(),
                       isSelected.data());

      // We need to store all the calculated uv values as we will need them when
      // apply the weights/corrections inline
      std::copy_n(scratch_image_weights_.data(), selectedBand.ChannelCount(),
                  uv_buffer_index);

      // Prime some of the required correction caches so that the correct data
      // is already available when we calculate/apply them inline
      PreCacheCorrections(*meta_data, msData.antenna_names, selectedBand);
      already_counted_weights.push_back(false);

      // We don't currently store the imaging weights
      // For testing we don't need this
      // Final implementation will need this
      // if (StoreImagingWeights())
      // ms_reader.WriteImagingWeights(scratch_image_weights_.data());

      weight_buffer_index += data_size;
      vis_buffer_index += data_size;
      uv_buffer_index += selectedBand.ChannelCount();
      uvw_buffer_index += 3;
      ++meta_data;

      ++nRows;
      msReader->NextInputRow();
    }

    std::function<std::complex<float>&(size_t index)> getVisibilityWithWeight;
    if (n_vis_polarizations == 1) {
      getVisibilityWithWeight = [&](size_t index) -> std::complex<float>& {
        size_t row = index / data_size;
        size_t channel = index % data_size;
        MSProvider::MetaData& meta_data = metadata_buffer[row];

        return InlineApplyWeightsAndCorrections<1>(
            msData.antenna_names, selectedBand, uv_buffer.data(),
            weight_buffer.data(), visibility_buffer.data(), index, row, channel,
            meta_data);
      };
    }

    auto ms = cmav<std::complex<float>, 2, cvirtmembuf<std::complex<float>>>(
        {nRows, selectedBand.ChannelCount()},
        cvirtmembuf<std::complex<float>>(getVisibilityWithWeight));

    Logger::Info << "Gridding " << nRows << " rows...\n";
    gridder_->AddInversionData(nRows, selectedBand.ChannelCount(),
                               uvw_buffer.data(), frequencies.data(), ms);

    // Clear 'per chunk' data that we have cached
    time_offsets_.clear();
    already_counted_weights.clear();

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
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getActualTrimmedSize(trimmedWidth, trimmedHeight);

  gridder_ = MakeGridder(trimmedWidth, trimmedHeight);
  gridder_->InitializeInversion();

  resetVisibilityCounters();

  for (MSData& msData : msDataVector) {
    gridMeasurementSet(msData);
  }

  gridder_->FinalizeImage(1.0 / totalWeight());

  Logger::Info << "Gridded visibility count: "
               << double(GriddedVisibilityCount());
  if (Weighting().IsNatural())
    Logger::Info << ", effective count after weighting: "
                 << EffectiveGriddedVisibilityCount();
  Logger::Info << '\n';

  _image = Image(ActualInversionWidth(), ActualInversionHeight());
  {
    std::vector<float> imageFloat = gridder_->RealImage();
    for (size_t i = 0; i < imageFloat.size(); ++i) _image[i] = imageFloat[i];
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    // Interpolate the image
    // The input is of size ActualInversionWidth() x ActualInversionHeight()
    schaapcommon::fft::Resampler resampler(
        ActualInversionWidth(), ActualInversionHeight(), ImageWidth(),
        ImageHeight(), resources_.NCpus());

    Image resized(ImageWidth(), ImageHeight());
    resampler.Resample(_image.Data(), resized.Data());
    _image = std::move(resized);
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight()
                  << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';

    _image = _image.Trim(TrimWidth(), TrimHeight());
  }
}

void WGriddingMSGridder::Predict(std::vector<Image>&& images) {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

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
