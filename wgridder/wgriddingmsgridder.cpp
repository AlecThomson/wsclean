#include "wgriddingmsgridder.h"

#include "wgriddinggridder_simple.h"

#include "../gridding/msgriddermanager.h"
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

namespace wsclean {

WGriddingMSGridder::WGriddingMSGridder(
    const Settings& settings, const Resources& resources,
    MsProviderCollection& ms_provider_collection, bool use_tuned_wgridder)
    : MsGridder(settings, ms_provider_collection),
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

size_t WGriddingMSGridder::CalculateConstantMemory() const {
  size_t constant_mem = gridder_->ConstantMemoryUsage();
  constant_mem += GetVisibilityModifier().GetCacheParmResponseSize();
  return constant_mem;
}

size_t WGriddingMSGridder::CalculateMaxRowsInMemory(
    int64_t available_memory, size_t constant_memory,
    size_t additional_per_row_consumption, size_t data_size) const {
  if (static_cast<int64_t>(constant_memory) >= available_memory) {
    // Assume that half the memory is necessary for the constant parts (like
    // image grid), and the other half remains available for the dynamic buffers
    constant_memory = available_memory / 2;
    Logger::Warn << "Not enough memory available for doing the gridding:\n"
                    "swapping might occur!\n";
  }

  const size_t per_visibility_ducc_overhead =
      gridder_->PerVisibilityMemoryUsage();
  const size_t per_row_visibility_memory =
      (per_visibility_ducc_overhead + sizeof(std::complex<float>)) * data_size;
  const size_t per_row_uvw_memory = sizeof(double) * 3;
  const uint64_t memory_per_row =
      additional_per_row_consumption  // external overheads
      + per_row_visibility_memory     // visibilities
      + per_row_uvw_memory;           // uvw
  const uint64_t memory_for_buffers = available_memory - constant_memory;
  const size_t max_n_rows =
      std::max(memory_for_buffers / memory_per_row, uint64_t(100));
  if (max_n_rows < 1000) {
    Logger::Warn << "Less than 1000 data rows fit in memory: this probably "
                    "means performance is going to be very poor!\n";
  }

  return max_n_rows;
}

size_t WGriddingMSGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();
  const aocommon::BandData selected_band(ms_data.SelectedBand());

  const size_t data_size = selected_band.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(data_size);
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<bool> selection_buffer(selected_band.ChannelCount(), true);

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  size_t max_rows_per_chunk = CalculateMaxRowsInMemory(
      resources_.Memory(), CalculateConstantMemory(), 0, data_size);

  aocommon::UVector<std::complex<float>> visibility_buffer(
      max_rows_per_chunk * selected_band.ChannelCount());
  aocommon::UVector<double> uvw_buffer(max_rows_per_chunk * 3);

  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  aocommon::UVector<std::complex<float>> row_visibilities(data_size);
  InversionRow row_data;
  row_data.data = row_visibilities.data();

  // Iterate over chunks until all data has been gridded
  size_t n_total_rows_read = 0;
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Debug << "Max " << max_rows_per_chunk << " rows fit in memory.\n";
    Logger::Info << "Loading data in memory...\n";

    size_t n_chunk_rows_read = 0;

    // Read / fill the chunk
    while (ms_reader->CurrentRowAvailable() &&
           n_chunk_rows_read < max_rows_per_chunk) {
      MSProvider::MetaData metadata;
      ms_reader->ReadMeta(metadata);
      row_data.uvw[0] = metadata.uInM;
      row_data.uvw[1] = metadata.vInM;
      row_data.uvw[2] = metadata.wInM;

      GetCollapsedVisibilities(*ms_reader, ms_data.antenna_names.size(),
                               row_data, selected_band, weight_buffer.data(),
                               model_buffer.data(), selection_buffer.data(),
                               metadata);

      std::copy_n(
          row_data.data, selected_band.ChannelCount(),
          &visibility_buffer[n_chunk_rows_read * selected_band.ChannelCount()]);
      std::copy_n(row_data.uvw, 3, &uvw_buffer[n_chunk_rows_read * 3]);

      ++n_chunk_rows_read;
      ms_reader->NextInputRow();
    }

    Logger::Info << "Gridding " << n_chunk_rows_read << " rows...\n";
    gridder_->AddInversionData(n_chunk_rows_read, selected_band.ChannelCount(),
                               uvw_buffer.data(), frequencies.data(),
                               visibility_buffer.data());

    n_total_rows_read += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_read;
}

size_t WGriddingMSGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  ms_data.ms_provider->ReopenRW();
  const aocommon::BandData selected_band(ms_data.SelectedBand());

  size_t n_total_rows_read = 0;

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  size_t max_rows_per_chunk =
      CalculateMaxRowsInMemory(resources_.Memory(), CalculateConstantMemory(),
                               0, selected_band.ChannelCount());

  aocommon::UVector<double> uvw_buffer(max_rows_per_chunk * 3);
  // Iterate over chunks until all data has been gridded
  ms_data.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    size_t n_chunk_rows_read = 0;

    // Read / fill the chunk
    Logger::Info << "Loading metadata...\n";
    // Read from metadata buffer
    std::vector<MSProvider::MetaData> metadata_buffer;
    while (ms_reader->CurrentRowAvailable() &&
           n_chunk_rows_read < max_rows_per_chunk) {
      MSProvider::MetaData metadata;
      ReadPredictMetaData(metadata);
      uvw_buffer[n_chunk_rows_read * 3] = metadata.uInM;
      uvw_buffer[n_chunk_rows_read * 3 + 1] = metadata.vInM;
      uvw_buffer[n_chunk_rows_read * 3 + 2] = metadata.wInM;
      metadata_buffer.emplace_back(std::move(metadata));
      n_chunk_rows_read++;

      ms_reader->NextInputRow();
    }

    Logger::Info << "Predicting " << n_chunk_rows_read << " rows...\n";
    aocommon::UVector<std::complex<float>> visibility_buffer(
        max_rows_per_chunk * selected_band.ChannelCount());
    gridder_->PredictVisibilities(
        n_chunk_rows_read, selected_band.ChannelCount(), uvw_buffer.data(),
        frequencies.data(), visibility_buffer.data());

    Logger::Info << "Writing...\n";
    for (size_t row = 0; row != n_chunk_rows_read; ++row) {
      WriteCollapsedVisibilities(
          *ms_data.ms_provider, ms_data.antenna_names.size(), selected_band,
          &visibility_buffer[row * selected_band.ChannelCount()],
          metadata_buffer[row]);
    }
    n_total_rows_read += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_read;
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

void WGriddingMSGridder::StartInversion() {
  size_t trimmed_width;
  size_t trimmed_height;
  getActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);
  gridder_->InitializeInversion();

  ResetVisibilityCounters();
}

void WGriddingMSGridder::FinishInversion() {
  gridder_->FinalizeImage(1.0 / ImageWeight());

  std::string log_message =
      "Gridded visibility count: " + std::to_string(GriddedVisibilityCount());
  if (Weighting().IsNatural()) {
    log_message += ", effective count after weighting: " +
                   std::to_string(EffectiveGriddedVisibilityCount());
  }
  Logger::Info << log_message + '\n';

  _image = Image(ActualInversionWidth(), ActualInversionHeight());
  {
    std::vector<float> image_float = gridder_->RealImage();
    for (size_t i = 0; i < image_float.size(); ++i) _image[i] = image_float[i];
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

void WGriddingMSGridder::StartPredict(std::vector<Image>&& images) {
  size_t trimmed_width;
  size_t trimmed_height;
  getActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmed_image(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    Image::Untrim(untrimmed_image.Data(), ImageWidth(), ImageHeight(),
                  images[0].Data(), TrimWidth(), TrimHeight());
    images[0] = std::move(untrimmed_image);
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    Image resampled_image(ImageWidth(), ImageHeight());
    schaapcommon::fft::Resampler resampler(
        ImageWidth(), ImageHeight(), ActualInversionWidth(),
        ActualInversionHeight(), resources_.NCpus());

    resampler.Resample(images[0].Data(), resampled_image.Data());
    images[0] = std::move(resampled_image);
  }

  gridder_->InitializePrediction(images[0].Data());
  images[0].Reset();
}

void WGriddingMSGridder::FinishPredict() {}

}  // namespace wsclean
