#include "msgridderbase.h"

#include "../math/calculatefftsize.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include "../structures/imageweights.h"

#include <aocommon/logger.h>
#include <aocommon/units/angle.h>

#include <limits>
#include <aocommon/multibanddata.h>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include <casacore/casa/Arrays/Cube.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/TableRecord.h>

using aocommon::Logger;
using schaapcommon::h5parm::JonesParameters;

namespace {

/**
 * @brief Select unique times from a given MSProvider
 */
std::vector<double> SelectUniqueTimes(MSProvider& ms_provider) {
  std::unique_ptr<MSReader> msReader = ms_provider.MakeReader();
  std::vector<double> msTimes;
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);
    // Assumes that the time instants in the MS are in ascending order.
    // In case this is violated, the returned vector will contain redundant
    // entries.
    if (msTimes.empty() || metaData.time != msTimes.back()) {
      msTimes.push_back(metaData.time);
    }
    msReader->NextInputRow();
  }
  return msTimes;
}
}  // namespace

// Defined out of class to allow the class the be used in a std::unique_ptr.
MSGridderBase::~MSGridderBase() = default;

MSGridderBase::MSGridderBase(const Settings& settings)
    : settings_(settings),
      w_grid_size_(settings.nWLayers),
      data_column_name_(settings.dataColumnName),
      small_inversion_(settings.smallInversion),
      w_limit_(settings.wLimit / 100.0),
      weighting_(settings.weightMode),
      visibility_weighting_mode_(settings.visibilityWeightingMode) {
#ifdef HAVE_EVERYBEAM
  visibility_modifier_.SetBeamInfo(settings.beamMode,
                                   settings.beamNormalisationMode);
#endif
}

std::vector<std::string> MSGridderBase::getAntennaNames(
    const casacore::MSAntenna& msAntenna) {
  const casacore::ScalarColumn<casacore::String> antennaNameColumn(
      msAntenna, msAntenna.columnName(casacore::MSAntenna::NAME));

  std::vector<std::string> antenna_names;
  antenna_names.reserve(antennaNameColumn.nrow());
  for (size_t i = 0; i < antennaNameColumn.nrow(); ++i) {
    antenna_names.push_back(antennaNameColumn(i));
  }
  return antenna_names;
}

void MSGridderBase::initializePointResponse(
    const MSGridderBase::MSData& msData) {
#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    const std::string element_response_string =
        !settings_.beamModel.empty() ? settings_.beamModel : "DEFAULT";
    visibility_modifier_.InitializePointResponse(
        msData.ms_provider->MS(), settings_.facetBeamUpdateTime,
        element_response_string, msData.bandData.ChannelCount(),
        settings_.dataColumnName, settings_.mwaPath);
  } else {
    visibility_modifier_.SetNoPointResponse();
  }
#else
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    throw std::runtime_error(
        "-apply-facet-beam or -grid-with-beam was set, but wsclean was not "
        "compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif
}

void MSGridderBase::StartMeasurementSet(const MSGridderBase::MSData& msData,
                                        bool isPredict) {
  initializePointResponse(msData);
  ms_index_ = msData.msIndex;
  n_vis_polarizations_ = msData.ms_provider->NPolarizations();
  gain_mode_ = GetGainMode(Polarization(), n_vis_polarizations_);
  const size_t n_channels = msData.SelectedBand().ChannelCount();
  scratch_image_weights_.resize(n_channels);
  if (isPredict) {
    scratch_model_data_.resize(n_channels *
                               msData.ms_provider->NPolarizations());
    predict_reader_ = msData.ms_provider->MakeReader();
  }
}

void MSGridderBase::initializeBandData(const casacore::MeasurementSet& ms,
                                       MSGridderBase::MSData& msData) {
  msData.bandData = aocommon::MultiBandData(ms)[msData.dataDescId];
  if (Selection(msData.msIndex).HasChannelRange()) {
    msData.startChannel = Selection(msData.msIndex).ChannelRangeStart();
    msData.endChannel = Selection(msData.msIndex).ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << msData.startChannel << '-'
                  << msData.endChannel << '\n';
    if (msData.startChannel >= msData.bandData.ChannelCount() ||
        msData.endChannel > msData.bandData.ChannelCount() ||
        msData.startChannel == msData.endChannel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << msData.bandData.ChannelCount()
          << " channels, requested imaging range is " << msData.startChannel
          << " -- " << msData.endChannel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    msData.startChannel = 0;
    msData.endChannel = msData.bandData.ChannelCount();
  }
}

template <size_t NPolInMSProvider>
void MSGridderBase::calculateWLimits(MSGridderBase::MSData& msData) {
  Logger::Info << "Determining min and max w & theoretical beam size... ";
  Logger::Info.Flush();
  msData.maxW = 0.0;
  msData.maxWWithFlags = 0.0;
  msData.minW = 1e100;
  msData.maxBaselineUVW = 0.0;
  msData.maxBaselineInM = 0.0;
  const aocommon::BandData selectedBand = msData.SelectedBand();
  std::vector<float> weightArray(selectedBand.ChannelCount() *
                                 NPolInMSProvider);
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();
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
    const double halfWidth = 0.5 * ImageWidth();
    const double halfHeight = 0.5 * ImageHeight();
    if (wHi > msData.maxW || wLo < msData.minW ||
        baselineInM / selectedBand.SmallestWavelength() >
            msData.maxBaselineUVW) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();
      for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
        const double wavelength = selectedBand.ChannelWavelength(ch);
        double wInL = metaData.wInM / wavelength;
        msData.maxWWithFlags = std::max(msData.maxWWithFlags, fabs(wInL));
        if (*weightPtr != 0.0) {
          double uInL = metaData.uInM / wavelength,
                 vInL = metaData.vInM / wavelength,
                 x = uInL * PixelSizeX() * ImageWidth(),
                 y = vInL * PixelSizeY() * ImageHeight(),
                 imagingWeight = GetImageWeights()->GetWeight(uInL, vInL);
          if (imagingWeight != 0.0) {
            if (std::floor(x) > -halfWidth && std::ceil(x) < halfWidth &&
                std::floor(y) > -halfHeight && std::ceil(y) < halfHeight) {
              msData.maxW = std::max(msData.maxW, std::fabs(wInL));
              msData.minW = std::min(msData.minW, std::fabs(wInL));
              msData.maxBaselineUVW =
                  std::max(msData.maxBaselineUVW, baselineInM / wavelength);
              msData.maxBaselineInM =
                  std::max(msData.maxBaselineInM, baselineInM);
            }
          }
        }
        weightPtr += NPolInMSProvider;
      }
    }

    msReader->NextInputRow();
  }

  if (msData.minW == 1e100) {
    msData.minW = 0.0;
    msData.maxWWithFlags = 0.0;
    msData.maxW = 0.0;
  }

  Logger::Info << "DONE (w=[" << msData.minW << ":" << msData.maxW
               << "] lambdas, maxuvw=" << msData.maxBaselineUVW << " lambda)\n";
  if (msData.maxWWithFlags != msData.maxW) {
    Logger::Debug << "Discarded data has higher w value of "
                  << msData.maxWWithFlags << " lambda.\n";
  }

  if (lastTime == firstTime || nTimesteps < 2)
    msData.integrationTime = 1;
  else
    msData.integrationTime = (lastTime - firstTime) / (nTimesteps - 1);
}

template void MSGridderBase::calculateWLimits<1>(MSGridderBase::MSData& msData);
template void MSGridderBase::calculateWLimits<2>(MSGridderBase::MSData& msData);
template void MSGridderBase::calculateWLimits<4>(MSGridderBase::MSData& msData);

void MSGridderBase::initializeMSDataVector(
    std::vector<MSGridderBase::MSData>& msDataVector) {
  if (MeasurementSetCount() == 0)
    throw std::runtime_error(
        "Something is wrong during inversion: no measurement sets given to "
        "inversion algorithm");
  msDataVector = std::vector<MSGridderBase::MSData>(MeasurementSetCount());

  resetMetaData();
  // FIXME: migrate data members to GriddingResult
  meta_data_cache_->correction_sum = 0.0;
  meta_data_cache_->h5_correction_sum = 0.0;

  bool hasCache = !meta_data_cache_->msDataVector.empty();
  if (!hasCache) meta_data_cache_->msDataVector.resize(MeasurementSetCount());

  if (!settings_.facetSolutionFiles.empty()) {
    visibility_modifier_.ResetCache(MeasurementSetCount(),
                                    settings_.facetSolutionFiles,
                                    settings_.facetSolutionTables);
  }

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    msDataVector[i].msIndex = i;
    initializeMeasurementSet(msDataVector[i], meta_data_cache_->msDataVector[i],
                             hasCache);
  }
  calculateOverallMetaData(msDataVector);
}

void MSGridderBase::initializeMeasurementSet(MSGridderBase::MSData& msData,
                                             MetaDataCache::Entry& cacheEntry,
                                             bool isCacheInitialized) {
  MSProvider& ms_provider = MeasurementSet(msData.msIndex);
  msData.ms_provider = &ms_provider;
  SynchronizedMS ms(ms_provider.MS());
  if (ms->nrow() == 0) throw std::runtime_error("Table has no rows (no data)");

  msData.antenna_names = getAntennaNames(ms->antenna());
  msData.dataDescId = ms_provider.DataDescId();

  initializeBandData(*ms, msData);

  if (HasDenormalPhaseCentre())
    Logger::Debug << "Set has denormal phase centre: dl=" << l_shift_
                  << ", dm=" << m_shift_ << '\n';

  calculateMSLimits(msData.SelectedBand(), ms_provider.StartTime());

  if (isCacheInitialized) {
    msData.maxW = cacheEntry.max_w;
    msData.maxWWithFlags = cacheEntry.max_w_with_flags;
    msData.minW = cacheEntry.min_w;
    msData.maxBaselineUVW = cacheEntry.max_baseline_uvw;
    msData.maxBaselineInM = cacheEntry.max_baseline_in_m;
    msData.integrationTime = cacheEntry.integration_time;
  } else {
    if (ms_provider.NPolarizations() == 4)
      calculateWLimits<4>(msData);
    else if (ms_provider.NPolarizations() == 2)
      calculateWLimits<2>(msData);
    else
      calculateWLimits<1>(msData);
    cacheEntry.max_w = msData.maxW;
    cacheEntry.max_w_with_flags = msData.maxWWithFlags;
    cacheEntry.min_w = msData.minW;
    cacheEntry.max_baseline_uvw = msData.maxBaselineUVW;
    cacheEntry.max_baseline_in_m = msData.maxBaselineInM;
    cacheEntry.integration_time = msData.integrationTime;
  }

  if (!settings_.facetSolutionFiles.empty()) {
    visibility_modifier_.SetMSTimes(msData.msIndex,
                                    SelectUniqueTimes(*msData.ms_provider));
  }
}

void MSGridderBase::calculateOverallMetaData(
    const std::vector<MSData>& msDataVector) {
  max_w_ = 0.0;
  min_w_ = std::numeric_limits<double>::max();
  double maxBaseline = 0.0;

  for (const MSData& msData : msDataVector) {
    maxBaseline = std::max(maxBaseline, msData.maxBaselineUVW);
    max_w_ = std::max(max_w_, msData.maxW);
    min_w_ = std::min(min_w_, msData.minW);
  }
  if (min_w_ > max_w_) {
    min_w_ = max_w_;
    Logger::Error
        << "*** Error! ***\n"
           "*** Calculating maximum and minimum w values failed! Make sure the "
           "data selection and scale settings are correct!\n"
           "***\n";
  }

  theoretical_beam_size_ = 1.0 / maxBaseline;
  if (IsFirstTask()) {
    Logger::Info << "Theoretic beam = "
                 << aocommon::units::Angle::ToNiceString(theoretical_beam_size_)
                 << "\n";
  }
  if (w_limit_ != 0.0) {
    max_w_ *= (1.0 - w_limit_);
    if (max_w_ < min_w_) max_w_ = min_w_;
  }

  if (!HasTrimSize()) SetTrimSize(ImageWidth(), ImageHeight());

  actual_inversion_width_ = ImageWidth();
  actual_inversion_height_ = ImageHeight();
  actual_pixel_size_x_ = PixelSizeX();
  actual_pixel_size_y_ = PixelSizeY();

  if (SmallInversion()) {
    size_t optWidth, optHeight, minWidth, minHeight;
    CalculateFFTSize(actual_inversion_width_, actual_pixel_size_x_,
                     theoretical_beam_size_, minWidth, optWidth);
    CalculateFFTSize(actual_inversion_height_, actual_pixel_size_y_,
                     theoretical_beam_size_, minHeight, optHeight);
    if (optWidth < actual_inversion_width_ ||
        optHeight < actual_inversion_height_) {
      const size_t newWidth =
          std::max(std::min(optWidth, actual_inversion_width_), size_t(32));
      const size_t newHeight =
          std::max(std::min(optHeight, actual_inversion_height_), size_t(32));
      if (IsFirstTask()) {
        Logger::Info << "Minimal inversion size: " << minWidth << " x "
                     << minHeight << ", using optimal: " << newWidth << " x "
                     << newHeight << "\n";
      }
      actual_pixel_size_x_ =
          (double(actual_inversion_width_) * actual_pixel_size_x_) /
          double(newWidth);
      actual_pixel_size_y_ =
          (double(actual_inversion_height_) * actual_pixel_size_y_) /
          double(newHeight);
      actual_inversion_width_ = newWidth;
      actual_inversion_height_ = newHeight;
    } else {
      if (IsFirstTask()) {
        Logger::Info
            << "Small inversion enabled, but inversion resolution already "
               "smaller than beam size: not using optimization.\n";
      }
    }
  }

  // Always call getSuggestedWGridSize in the first iteration, since it then
  // logs the suggested wgrid size.
  const size_t suggestedGridSize =
      (IsFirstTask() || !hasWGridSize()) ? getSuggestedWGridSize() : 0;
  actual_w_grid_size_ = hasWGridSize() ? w_grid_size_ : suggestedGridSize;
}

void MSGridderBase::ReadPredictMetaData(MSProvider::MetaData& metaData) {
  predict_reader_->ReadMeta(metaData);
  predict_reader_->NextInputRow();
}

template <size_t PolarizationCount, GainMode GainEntry>
void MSGridderBase::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData) {
  assert(GetPsfMode() == PsfMode::kNone);  // The PSF is never predicted.

  if (visibility_modifier_.HasH5Parm()) {
    assert(!settings_.facetRegionFilename.empty());
    visibility_modifier_.CacheParmResponse(metaData.time, antenna_names,
                                           curBand, ms_index_);

    visibility_modifier_.ApplyParmResponse<PolarizationCount, GainEntry>(
        buffer, ms_index_, curBand.ChannelCount(), antenna_names.size(),
        metaData.antenna1, metaData.antenna2);
  }

#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam) {
    visibility_modifier_.CacheBeamResponse(metaData.time, metaData.fieldId,
                                           curBand);

    visibility_modifier_.ApplyBeamResponse<PolarizationCount, GainEntry>(
        buffer, curBand.ChannelCount(), metaData.antenna1, metaData.antenna2);
  }
#endif

  {
    const size_t lock_index =
        facet_group_index_ * MeasurementSetCount() + ms_index_;
    std::unique_ptr<GriddingTaskManager::WriterLock> lock =
        writer_lock_manager_->GetLock(lock_index);
    ms_provider.WriteModel(buffer, IsFacet());
  }
  ms_provider.NextOutputRow();
}

template void MSGridderBase::WriteInstrumentalVisibilities<1, GainMode::kXX>(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void MSGridderBase::WriteInstrumentalVisibilities<1, GainMode::kYY>(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void
MSGridderBase::WriteInstrumentalVisibilities<1, GainMode::kDiagonal>(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void
MSGridderBase::WriteInstrumentalVisibilities<2, GainMode::kDiagonal>(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void MSGridderBase::WriteInstrumentalVisibilities<4, GainMode::kFull>(
    MSProvider& ms_provider, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template <size_t PolarizationCount>
void MSGridderBase::CalculateWeights(InversionRow& rowData,
                                     const aocommon::BandData& curBand,
                                     float* weightBuffer,
                                     std::complex<float>* modelBuffer,
                                     const bool* isSelected) {
  const std::size_t dataSize = curBand.ChannelCount() * PolarizationCount;
  if (GetPsfMode() != PsfMode::kNone) {
    // Visibilities for a point source at the phase centre are all ones
    std::fill_n(rowData.data, dataSize, 1.0);
    double dl = 0.0;
    double dm = 0.0;
    if (GetPsfMode() == PsfMode::kSingle) {
      // The point source is shifted to the centre of the main image
      dl = MainImageDL();
      dm = MainImageDM();
    } else {  // GetPsfMode() == PsfMode::kDirectionDependent
      // The point source is shifted to the centre of the current DdPsf position
      dl = LShift();
      dm = MShift();
    }
    if (dl != 0.0 || dm != 0.0) {
      const double dn = std::sqrt(1.0 - dl * dl - dm * dm) - 1.0;
      const double shiftFactor =
          2.0 * M_PI *
          (rowData.uvw[0] * dl + rowData.uvw[1] * dm + rowData.uvw[2] * dn);
      rotateVisibilities<PolarizationCount>(curBand, shiftFactor, rowData.data);
    }
  }

  if (DoSubtractModel()) {
    std::complex<float>* modelIter = modelBuffer;
    for (std::complex<float>* iter = rowData.data;
         iter != rowData.data + dataSize; ++iter) {
      *iter -= *modelIter;
      modelIter++;
    }
  }

  // Any visibilities that are not gridded in this pass
  // should not contribute to the weight sum, so set these
  // to have zero weight.
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    for (size_t p = 0; p != PolarizationCount; ++p) {
      if (!isSelected[ch]) weightBuffer[ch * PolarizationCount + p] = 0.0;
    }
  }

  switch (GetVisibilityWeightingMode()) {
    case VisibilityWeightingMode::NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case VisibilityWeightingMode::SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t i = 0; i != dataSize; ++i) weightBuffer[i] *= weightBuffer[i];
      break;
    case VisibilityWeightingMode::UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t i = 0; i != dataSize; ++i) {
        if (weightBuffer[i] != 0.0) weightBuffer[i] = 1.0f;
      }
      break;
  }

  // Precompute imaging weights
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    const double u = rowData.uvw[0] / curBand.ChannelWavelength(ch);
    const double v = rowData.uvw[1] / curBand.ChannelWavelength(ch);
    scratch_image_weights_[ch] = GetImageWeights()->GetWeight(u, v);
  }
}

template <size_t PolarizationCount, GainMode GainEntry>
void MSGridderBase::ApplyWeightsAndCorrections(
    const std::vector<std::string>& antenna_names, InversionRow& rowData,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData) {
  if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    const bool apply_forward = GetPsfMode() == PsfMode::kDirectionDependent;
    if (apply_beam && visibility_modifier_.HasH5Parm()) {
#ifdef HAVE_EVERYBEAM
      // Load and apply (in conjugate) both the beam and the h5parm solutions
      visibility_modifier_.CacheBeamResponse(metaData.time, metaData.fieldId,
                                             curBand);
      visibility_modifier_.CacheParmResponse(metaData.time, antenna_names,
                                             curBand, ms_index_);
      const VisibilityModifier::DualResult result =
          visibility_modifier_
              .ApplyConjugatedDual<PolarizationCount, GainEntry>(
                  rowData.data, weightBuffer, scratch_image_weights_.data(),
                  curBand.ChannelCount(), antenna_names.size(),
                  metaData.antenna1, metaData.antenna2, ms_index_,
                  apply_forward);
      meta_data_cache_->h5_correction_sum += result.h5Sum;
      meta_data_cache_->correction_sum += result.correctionSum;
    } else if (apply_beam) {
      // Load and apply only the conjugate beam
      visibility_modifier_.CacheBeamResponse(metaData.time, metaData.fieldId,
                                             curBand);
      meta_data_cache_->correction_sum +=
          visibility_modifier_
              .ApplyConjugatedBeamResponse<PolarizationCount, GainEntry>(
                  rowData.data, weightBuffer, scratch_image_weights_.data(),
                  curBand.ChannelCount(), metaData.antenna1, metaData.antenna2,
                  apply_forward);

#endif  // HAVE_EVERYBEAM
    } else if (visibility_modifier_.HasH5Parm()) {
      visibility_modifier_.CacheParmResponse(metaData.time, antenna_names,
                                             curBand, ms_index_);

      meta_data_cache_->correction_sum +=
          visibility_modifier_
              .ApplyConjugatedParmResponse<PolarizationCount, GainEntry>(
                  rowData.data, weightBuffer, scratch_image_weights_.data(),
                  ms_index_, curBand.ChannelCount(), antenna_names.size(),
                  metaData.antenna1, metaData.antenna2, apply_forward);
    }
  }

  // Apply visibility and imaging weights
  std::complex<float>* dataIter = rowData.data;
  float* weightIter = weightBuffer;
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    for (size_t p = 0; p != PolarizationCount; ++p) {
      const double cumWeight = *weightIter * scratch_image_weights_[ch];
      if (p == 0 && cumWeight != 0.0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        visibility_weight_sum_ += *weightIter;
        max_gridded_weight_ = std::max(cumWeight, max_gridded_weight_);
        ++gridded_visibility_count_;
      }
      // Total weight includes imaging weights
      total_weight_ += cumWeight;
      *weightIter = cumWeight;
      *dataIter *= cumWeight;
      ++dataIter;
      ++weightIter;
    }
  }
}

template void MSGridderBase::CalculateWeights<1>(
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);
template void MSGridderBase::CalculateWeights<2>(
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);
template void MSGridderBase::CalculateWeights<3>(
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);
template void MSGridderBase::CalculateWeights<4>(
    InversionRow& newItem, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::ApplyWeightsAndCorrections<1, GainMode::kXX>(
    const std::vector<std::string>& antenna_names, InversionRow& newItem,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData);

template void MSGridderBase::ApplyWeightsAndCorrections<1, GainMode::kYY>(
    const std::vector<std::string>& antenna_names, InversionRow& newItem,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData);

template void MSGridderBase::ApplyWeightsAndCorrections<1, GainMode::kDiagonal>(
    const std::vector<std::string>& antenna_names, InversionRow& newItem,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData);

template void MSGridderBase::ApplyWeightsAndCorrections<2, GainMode::kDiagonal>(
    const std::vector<std::string>& antenna_names, InversionRow& newItem,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData);

template void MSGridderBase::ApplyWeightsAndCorrections<4, GainMode::kFull>(
    const std::vector<std::string>& antenna_names, InversionRow& newItem,
    const aocommon::BandData& curBand, float* weightBuffer,
    const MSProvider::MetaData& metaData);

template <size_t PolarizationCount>
void MSGridderBase::rotateVisibilities(const aocommon::BandData& bandData,
                                       double shiftFactor,
                                       std::complex<float>* dataIter) {
  for (size_t ch = 0; ch != bandData.ChannelCount(); ++ch) {
    const double wShiftRad = shiftFactor / bandData.ChannelWavelength(ch);
    const std::complex<float> phasor(std::cos(wShiftRad), std::sin(wShiftRad));
    for (size_t p = 0; p != PolarizationCount; ++p) {
      *dataIter *= phasor;
      ++dataIter;
    }
  }
}

template void MSGridderBase::rotateVisibilities<1>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<2>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<4>(
    const aocommon::BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
