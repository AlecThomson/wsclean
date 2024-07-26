#include "msprovidercollection.h"

#include <vector>

#include <aocommon/logger.h>

#include "msgridderbase.h"
#include "msgriddermanager.h"

using aocommon::Logger;

namespace wsclean {

namespace {
template <size_t NPolInMSProvider>
inline void CalculateMsLimits(MsProviderCollection::MsData& ms_data, double u,
                              double v, double w, double baseline_in_meters,
                              double wavelength, double pixel_size_x,
                              double pixel_size_y, size_t image_width,
                              size_t image_height,
                              const ImageWeights* image_weights) {
  const double half_width = 0.5 * image_width;
  const double half_height = 0.5 * image_height;
  const double x = u * pixel_size_x * image_width;
  const double y = u * pixel_size_y * image_height;
  const double imaging_weight = image_weights->GetWeight(u, v);
  if (imaging_weight != 0.0) {
    if (std::floor(x) > -half_width && std::ceil(x) < half_width &&
        std::floor(y) > -half_height && std::ceil(y) < half_height) {
      ms_data.max_w = std::max(ms_data.max_w, std::fabs(w));
      ms_data.min_w = std::min(ms_data.min_w, std::fabs(w));
      ms_data.max_baseline_uvw =
          std::max(ms_data.max_baseline_uvw, baseline_in_meters / wavelength);
      ms_data.max_baseline_meters =
          std::max(ms_data.max_baseline_meters, baseline_in_meters);
    }
  }
}
}  // namespace

void MsProviderCollection::InitializeMSDataVector(
    const std::vector<MSGridderBase*>& gridders, double w_limit) {
  assert(Count() != 0);

  bool hasCache = false;
  for (MSGridderBase* facet_gridder : gridders) {
    hasCache = facet_gridder->HasMetaDataCache();
    if (!hasCache) facet_gridder->AllocateMetaDataCache(Count());

    facet_gridder->ResetVisibilityModifierCache();
  }

  ms_limits_.max_baseline = 0.0;
  ms_limits_.max_w = 0.0;
  ms_limits_.min_w = std::numeric_limits<double>::max();
  for (size_t i = 0; i != Count(); ++i) {
    MsData& ms_data = ms_data_vector_[i];
    ms_data.internal_ms_index = i;
    ms_data.original_ms_index = Index(i);
    InitializeMeasurementSet(ms_data, gridders, hasCache);

    ms_limits_.Calculate(ms_data.SelectedBand(),
                         ms_data.ms_provider->StartTime());
    ms_limits_.max_baseline =
        std::max(ms_limits_.max_baseline, ms_data.max_baseline_uvw);
    ms_limits_.max_w = std::max(ms_limits_.max_w, ms_data.max_w);
    ms_limits_.min_w = std::min(ms_limits_.min_w, ms_data.min_w);
  }

  ms_limits_.Validate();

  if (w_limit != 0.0) {
    ms_limits_.max_w *= (1.0 - w_limit);
    if (ms_limits_.max_w < ms_limits_.min_w)
      ms_limits_.max_w = ms_limits_.min_w;
  }

  for (MSGridderBase* facet_gridder : gridders) {
    facet_gridder->SetMaxW(ms_limits_.max_w);
    facet_gridder->SetMinW(ms_limits_.min_w);
    facet_gridder->SetMaxBaseline(ms_limits_.max_baseline);
  }
}

void MsProviderCollection::InitializeMS() { ms_data_vector_.resize(Count()); }

std::vector<std::string> MsProviderCollection::GetAntennaNames(
    const casacore::MSAntenna& antenna) {
  const casacore::ScalarColumn<casacore::String> antennaNameColumn(
      antenna, antenna.columnName(casacore::MSAntenna::NAME));

  std::vector<std::string> antenna_names;
  antenna_names.reserve(antennaNameColumn.nrow());
  for (size_t i = 0; i < antennaNameColumn.nrow(); ++i) {
    antenna_names.push_back(antennaNameColumn(i));
  }
  return antenna_names;
}

void MsProviderCollection::MsData::InitializeBandData(
    const casacore::MeasurementSet& ms, const MSSelection& selection) {
  band_data = aocommon::MultiBandData(ms)[data_desc_id];
  if (selection.HasChannelRange()) {
    start_channel = selection.ChannelRangeStart();
    end_channel = selection.ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << start_channel << '-'
                  << end_channel << '\n';
    if (start_channel >= band_data.ChannelCount() ||
        end_channel > band_data.ChannelCount() ||
        start_channel == end_channel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << band_data.ChannelCount()
          << " channels, requested imaging range is " << start_channel << " -- "
          << end_channel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    start_channel = 0;
    end_channel = band_data.ChannelCount();
  }
}

void MsProviderCollection::InitializeMeasurementSet(
    MsData& ms_data, const std::vector<MSGridderBase*>& gridders,
    bool is_cached) {
  MSProvider& ms_provider = MeasurementSet(ms_data.internal_ms_index);
  ms_data.ms_provider = &ms_provider;

  {
    SynchronizedMS ms(ms_provider.MS());
    if (ms->nrow() == 0)
      throw std::runtime_error("Table has no rows (no data)");

    ms_data.antenna_names = GetAntennaNames(ms->antenna());
    ms_data.data_desc_id = ms_provider.DataDescId();

    ms_data.InitializeBandData(*ms, Selection(ms_data.internal_ms_index));
  }

  // wlimits will vary across facets in a facet group, however these limits are
  // "estimates".
  // They should not be too small but can be too large (though this can have
  // performance implications.
  // As a code simplification and performance improvement we select the smallest
  // width and height out of all facets in a facet group to calculate the
  // wlimits rather than doing it individually for each one.
  size_t min_image_width = gridders[0]->ImageWidth();
  size_t min_image_height = gridders[0]->ImageHeight();
  for (const MSGridderBase* gridder : gridders) {
    min_image_width = std::min(min_image_width, gridder->ImageWidth());
    min_image_height = std::min(min_image_height, gridder->ImageHeight());
  }

  MetaDataCache::Entry& cache_entry =
      gridders[0]->GetMetaDataCacheItem(ms_data.internal_ms_index);

  if (is_cached) {
    ms_data.max_w = cache_entry.max_w;
    ms_data.max_w_with_flags = cache_entry.max_w_with_flags;
    ms_data.min_w = cache_entry.min_w;
    ms_data.max_baseline_uvw = cache_entry.max_baseline_uvw;
    ms_data.max_baseline_meters = cache_entry.max_baseline_in_m;
    ms_data.integration_time = cache_entry.integration_time;
  } else {
    if (ms_provider.NPolarizations() == 4)
      CalculateMsLimits<4>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    else if (ms_provider.NPolarizations() == 2)
      CalculateMsLimits<2>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    else
      CalculateMsLimits<1>(ms_data, gridders[0]->PixelSizeX(),
                           gridders[0]->PixelSizeY(), min_image_width,
                           min_image_height, gridders[0]->GetImageWeights());
    cache_entry.max_w = ms_data.max_w;
    cache_entry.max_w_with_flags = ms_data.max_w_with_flags;
    cache_entry.min_w = ms_data.min_w;
    cache_entry.max_baseline_uvw = ms_data.max_baseline_uvw;
    cache_entry.max_baseline_in_m = ms_data.max_baseline_meters;
    cache_entry.integration_time = ms_data.integration_time;
  }

  for (MSGridderBase* gridder : gridders) {
    gridder->initializeVisibilityModifierTimes(ms_data);
  }
}

template <size_t NPolInMSProvider>
void MsProviderCollection::CalculateMsLimits(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height,
    const ImageWeights* image_weights) {
  Logger::Info << "Determining min and max w & theoretical beam size... ";
  Logger::Info.Flush();
  ms_data.max_w = 0.0;
  ms_data.max_w_with_flags = 0.0;
  ms_data.min_w = 1e100;
  ms_data.max_baseline_uvw = 0.0;
  ms_data.max_baseline_meters = 0.0;
  const aocommon::BandData selectedBand = ms_data.SelectedBand();
  std::vector<float> weightArray(selectedBand.ChannelCount() *
                                 NPolInMSProvider);
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  std::unique_ptr<MSReader> msReader = ms_data.ms_provider->MakeReader();
  const double smallestWavelength = selectedBand.SmallestWavelength();
  const double longestWavelength = selectedBand.LongestWavelength();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);

    if (curTimestep != metaData.time) {
      curTimestep = metaData.time;
      ++nTimesteps;
      if (firstTime == -1) firstTime = curTimestep;
      lastTime = curTimestep;
    }

    const double wHi = std::fabs(metaData.wInM / smallestWavelength);
    const double wLo = std::fabs(metaData.wInM / longestWavelength);
    const double baselineInM = std::sqrt(metaData.uInM * metaData.uInM +
                                         metaData.vInM * metaData.vInM +
                                         metaData.wInM * metaData.wInM);
    if (wHi > ms_data.max_w || wLo < ms_data.min_w ||
        baselineInM / selectedBand.SmallestWavelength() >
            ms_data.max_baseline_uvw) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();

      for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
        const double wavelength = selectedBand.ChannelWavelength(ch);
        double wInL = metaData.wInM / wavelength;
        ms_data.max_w_with_flags =
            std::max(ms_data.max_w_with_flags, fabs(wInL));
        if (*weightPtr != 0.0) {
          double uInL = metaData.uInM / wavelength;
          double vInL = metaData.vInM / wavelength;
          wsclean::CalculateMsLimits<NPolInMSProvider>(
              ms_data, uInL, vInL, wInL, baselineInM, wavelength, pixel_size_x,
              pixel_size_y, image_width, image_height, image_weights);
        }
        weightPtr += NPolInMSProvider;
      }
    }
    msReader->NextInputRow();
  }

  if (ms_data.min_w == 1e100) {
    ms_data.min_w = 0.0;
    ms_data.max_w_with_flags = 0.0;
    ms_data.max_w = 0.0;
  }

  Logger::Info << "DONE (w=[" << ms_data.min_w << ":" << ms_data.max_w
               << "] lambdas, maxuvw=" << ms_data.max_baseline_uvw
               << " lambda)\n";
  if (ms_data.max_w_with_flags != ms_data.max_w) {
    Logger::Debug << "Discarded data has higher w value of "
                  << ms_data.max_w_with_flags << " lambda.\n";
  }

  if (lastTime == firstTime || nTimesteps < 2)
    ms_data.integration_time = 1;
  else
    ms_data.integration_time = (lastTime - firstTime) / (nTimesteps - 1);
}

template void MsProviderCollection::CalculateMsLimits<1>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);
template void MsProviderCollection::CalculateMsLimits<2>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);
template void MsProviderCollection::CalculateMsLimits<4>(
    MsData& ms_data, double pixel_size_x, double pixel_size_y,
    size_t image_width, size_t image_height, const ImageWeights* image_weights);

}  // namespace wsclean
