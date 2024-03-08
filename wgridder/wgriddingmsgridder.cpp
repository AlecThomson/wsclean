
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
  const size_t dataSize = selectedBand.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> modelBuffer(dataSize);
  aocommon::UVector<float> weightBuffer(dataSize);
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  size_t totalNRows = 0;
  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  for (auto& gd : griddingData) {
    gd.nRows = 0;
    gd.lastChunk = false;
    gd.doneReading = false;
    gd.doneComputing = true;
    gd.uvwBuffer.resize(maxNRows * 3);
    gd.visBuffer.resize(maxNRows * selectedBand.ChannelCount());
  }

  MSProvider* msProvider = msData.ms_provider;
  std::vector<std::basic_string<char>>* antennaNames = &msData.antenna_names;

  std::complex<float>* f_modelBuffer = modelBuffer.data();
  float* f_weightBuffer = weightBuffer.data();
  bool* f_isSelected = isSelected.data();

  std::thread reader_thread([this, msProvider, selectedBand, f_modelBuffer,
                             f_weightBuffer, f_isSelected, maxNRows, dataSize,
                             antennaNames]() {
    this->read_fn(msProvider, selectedBand, f_modelBuffer, f_weightBuffer,
                  f_isSelected, maxNRows, dataSize, antennaNames);
  });
  // Iterate over chunks until all data has been gridded
  size_t chunk_idx = 0;
  while (true) {
    Logger::Debug << "Max " << maxNRows << " rows fit in memory.\n";
    Logger::Info << "Loading data in memory...\n";

    GriddingData& chunk = griddingData[chunk_idx % GRIDDING_PREFETCH_FACTOR];

    // Mutex wait for data
    if (!chunk.doneReading) {
      std::cout << "Computer is waiting for read." << std::endl;
    }
    while (!chunk.doneReading) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    Logger::Info << "Gridding " << chunk.nRows << " rows...\n";
    gridder_->AddInversionData(chunk.nRows, selectedBand.ChannelCount(),
                               chunk.uvwBuffer.data(), frequencies.data(),
                               chunk.visBuffer.data());

    totalNRows += chunk.nRows;
    bool done = chunk.lastChunk;

    chunk.doneReading = false;
    chunk.doneComputing = true;

    if (done) {
      break;
    }

    chunk_idx += 1;
  }  // end of chunk
  reader_thread.join();
  std::cout << "Fully done." << std::endl;

  msData.totalRowsProcessed += totalNRows;
}
void WGriddingMSGridder::read_fn(
    MSProvider* msProvider, const aocommon::BandData& selectedBand,
    std::complex<float>* modelBuffer, float* weightBuffer, bool* isSelected,
    size_t maxNRows, const size_t dataSize,
    std::vector<std::basic_string<char>>* antennaNames) {
  // Read / fill the chunk
  auto msReader = msProvider->MakeReader();
  aocommon::UVector<std::complex<float>> newItemData(dataSize);
  InversionRow newRowData;
  newRowData.data = newItemData.data();

  size_t chunk_idx = 0;
  while (true) {
    GriddingData& chunk = griddingData[chunk_idx % GRIDDING_PREFETCH_FACTOR];
    chunk.nRows = 0;

    // Mutex wait for data
    if (!chunk.doneComputing) {
    }

    while (!chunk.doneComputing) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    while (msReader->CurrentRowAvailable() && chunk.nRows < maxNRows) {
      MSProvider::MetaData metaData;
      msReader->ReadMeta(metaData);
      newRowData.uvw[0] = metaData.uInM;
      newRowData.uvw[1] = metaData.vInM;
      newRowData.uvw[2] = metaData.wInM;

      GetCollapsedVisibilities(*msReader, *antennaNames, newRowData,
                                     selectedBand, weightBuffer, modelBuffer,
                                     isSelected, metaData);

      std::copy_n(newRowData.data, selectedBand.ChannelCount(),
                  &chunk.visBuffer[chunk.nRows * selectedBand.ChannelCount()]);
      std::copy_n(newRowData.uvw, 3, &chunk.uvwBuffer[chunk.nRows * 3]);

      ++chunk.nRows;
      msReader->NextInputRow();
    }

    chunk.doneComputing = false;
    if (!msReader->CurrentRowAvailable()) {
      chunk.lastChunk = true;
      chunk.doneReading = true;
      break;
    }
    chunk.doneReading = true;

    chunk_idx += 1;
  }
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

    {
      // Lock over all writing operations
      // We do technically not need to lock here yet, but it saves
      const size_t lock_index =
          facet_group_index_ * MeasurementSetCount() + ms_index_;
      std::unique_ptr<GriddingTaskManager::WriterLock> lock =
          writer_lock_manager_->GetLock(lock_index);
      for (size_t row = 0; row != nRows; ++row) {
        FlushBufferVisibilities(*msData.ms_provider,
                                &visBuffer[row * selectedBand.ChannelCount()]);
      }
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
