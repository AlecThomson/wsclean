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
#include "wsmsgridder.h"

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"
#include "../main/settings.h"
#include "../structures/resources.h"
#include "../wgridder/wgriddingmsgridder.h"

using aocommon::Logger;

void MSManager::InitializeMSDataVector(
    const std::vector<MSGridderBase*>& gridders) {
  if (Count() == 0)
    throw std::runtime_error(
        "Something is wrong during inversion: no measurement sets given to "
        "inversion algorithm");

  ms_data_vector_.resize(Count());
  ms_facet_data_vector_.resize(gridders.size(),
                               std::vector<FacetData>(Count()));
  bool hasCache = false;
  for (auto& gridder : gridders) {
    hasCache = gridder->HasMetaDataCache();
    if (!hasCache) gridder->ReserveMetaDataCache(Count());

    gridder->ResetVisibilityModifierCache();
  }

  for (size_t i = 0; i != Count(); ++i) {
    ms_data_vector_[i].internal_ms_index = i;
    ms_data_vector_[i].original_ms_index = Index(i);
    InitializeMeasurementSet(ms_data_vector_[i], ms_facet_data_vector_,
                             gridders, hasCache);
  }
}

std::vector<std::string> MSManager::GetAntennaNames(
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

void MSManager::Data::InitializeBandData(const casacore::MeasurementSet& ms,
                                         const MSSelection& selection) {
  bandData = aocommon::MultiBandData(ms)[dataDescId];
  if (selection.HasChannelRange()) {
    startChannel = selection.ChannelRangeStart();
    endChannel = selection.ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << startChannel << '-' << endChannel
                  << '\n';
    if (startChannel >= bandData.ChannelCount() ||
        endChannel > bandData.ChannelCount() || startChannel == endChannel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << bandData.ChannelCount() << " channels, requested imaging range is "
          << startChannel << " -- " << endChannel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    startChannel = 0;
    endChannel = bandData.ChannelCount();
  }
}

void MSManager::InitializeMeasurementSet(
    Data& ms_data, std::vector<std::vector<FacetData>>& ms_facet_data_vector,
    const std::vector<MSGridderBase*>& gridders, bool is_cached) {
  MSProvider& ms_provider = MeasurementSet(ms_data.internal_ms_index);
  ms_data.ms_provider = &ms_provider;
  SynchronizedMS ms(ms_provider.MS());
  if (ms->nrow() == 0) throw std::runtime_error("Table has no rows (no data)");

  ms_data.antenna_names = GetAntennaNames(ms->antenna());
  ms_data.dataDescId = ms_provider.DataDescId();

  ms_data.InitializeBandData(*ms, Selection(ms_data.internal_ms_index));

  /*if (HasDenormalPhaseCentre())
    Logger::Debug << "Set has denormal phase centre: dl=" << l_shift_
                  << ", dm=" << m_shift_ << '\n';*/

  ms_limits_.Calculate(ms_data.SelectedBand(), ms_provider.StartTime());

  size_t gridder_index = 0;
  for (auto gridder : gridders) {
    MetaDataCache::Entry& cache_entry =
        gridder->GetMetaDataCacheItem(gridder_index);
    FacetData& facet_data =
        ms_facet_data_vector[gridder_index][ms_data.internal_ms_index];
    ++gridder_index;

    if (is_cached) {
      facet_data.maxW = cache_entry.max_w;
      facet_data.maxWWithFlags = cache_entry.max_w_with_flags;
      facet_data.minW = cache_entry.min_w;
      facet_data.maxBaselineUVW = cache_entry.max_baseline_uvw;
      facet_data.maxBaselineInM = cache_entry.max_baseline_in_m;
      facet_data.integrationTime = cache_entry.integration_time;
    } else {
      if (ms_provider.NPolarizations() == 4)
        CalculateWLimits<4>(facet_data, ms_data, gridder);
      else if (ms_provider.NPolarizations() == 2)
        CalculateWLimits<2>(facet_data, ms_data, gridder);
      else
        CalculateWLimits<1>(facet_data, ms_data, gridder);
      cache_entry.max_w = facet_data.maxW;
      cache_entry.max_w_with_flags = facet_data.maxWWithFlags;
      cache_entry.min_w = facet_data.minW;
      cache_entry.max_baseline_uvw = facet_data.maxBaselineUVW;
      cache_entry.max_baseline_in_m = facet_data.maxBaselineInM;
      cache_entry.integration_time = facet_data.integrationTime;
    }
    gridder->initializeVisibilityModifierTimes(ms_data);
  }
}

template <size_t NPolInMSProvider>
void MSManager::CalculateWLimits(FacetData& ms_facet_data, Data& ms_data,
                                 MSGridderBase* gridder) {
  Logger::Info << "Determining min and max w & theoretical beam size... ";
  Logger::Info.Flush();
  ms_facet_data.maxW = 0.0;
  ms_facet_data.maxWWithFlags = 0.0;
  ms_facet_data.minW = 1e100;
  ms_facet_data.maxBaselineUVW = 0.0;
  ms_facet_data.maxBaselineInM = 0.0;
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
    const double halfWidth = 0.5 * gridder->ImageWidth();
    const double halfHeight = 0.5 * gridder->ImageHeight();
    if (wHi > ms_facet_data.maxW || wLo < ms_facet_data.minW ||
        baselineInM / selectedBand.SmallestWavelength() >
            ms_facet_data.maxBaselineUVW) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();
      for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
        const double wavelength = selectedBand.ChannelWavelength(ch);
        double wInL = metaData.wInM / wavelength;
        ms_facet_data.maxWWithFlags =
            std::max(ms_facet_data.maxWWithFlags, fabs(wInL));
        if (*weightPtr != 0.0) {
          double uInL = metaData.uInM / wavelength,
                 vInL = metaData.vInM / wavelength,
                 x = uInL * gridder->PixelSizeX() * gridder->ImageWidth(),
                 y = vInL * gridder->PixelSizeY() * gridder->ImageHeight(),
                 imagingWeight =
                     gridder->GetImageWeights()->GetWeight(uInL, vInL);
          if (imagingWeight != 0.0) {
            if (std::floor(x) > -halfWidth && std::ceil(x) < halfWidth &&
                std::floor(y) > -halfHeight && std::ceil(y) < halfHeight) {
              ms_facet_data.maxW =
                  std::max(ms_facet_data.maxW, std::fabs(wInL));
              ms_facet_data.minW =
                  std::min(ms_facet_data.minW, std::fabs(wInL));
              ms_facet_data.maxBaselineUVW = std::max(
                  ms_facet_data.maxBaselineUVW, baselineInM / wavelength);
              ms_facet_data.maxBaselineInM =
                  std::max(ms_facet_data.maxBaselineInM, baselineInM);
            }
          }
        }
        weightPtr += NPolInMSProvider;
      }
    }

    msReader->NextInputRow();
  }

  if (ms_facet_data.minW == 1e100) {
    ms_facet_data.minW = 0.0;
    ms_facet_data.maxWWithFlags = 0.0;
    ms_facet_data.maxW = 0.0;
  }

  Logger::Info << "DONE (w=[" << ms_facet_data.minW << ":" << ms_facet_data.maxW
               << "] lambdas, maxuvw=" << ms_facet_data.maxBaselineUVW
               << " lambda)\n";
  if (ms_facet_data.maxWWithFlags != ms_facet_data.maxW) {
    Logger::Debug << "Discarded data has higher w value of "
                  << ms_facet_data.maxWWithFlags << " lambda.\n";
  }

  if (lastTime == firstTime || nTimesteps < 2)
    ms_facet_data.integrationTime = 1;
  else
    ms_facet_data.integrationTime = (lastTime - firstTime) / (nTimesteps - 1);
}

template void MSManager::CalculateWLimits<1>(FacetData& ms_facet_data,
                                             Data& ms_data,
                                             MSGridderBase* gridder);
template void MSManager::CalculateWLimits<2>(FacetData& ms_facet_data,
                                             Data& ms_data,
                                             MSGridderBase* gridder);
template void MSManager::CalculateWLimits<4>(FacetData& ms_facet_data,
                                             Data& ms_data,
                                             MSGridderBase* gridder);
