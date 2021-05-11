#include "msgridderbase.h"

#include "../io/logger.h"

#include "../math/calculatefftsize.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include "../structures/imageweights.h"

#include "../units/angle.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"
#include <limits>
#include <aocommon/matrix2x2.h>

#include <schaapcommon/h5parm/h5parm.h>
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

#include <atomic>

using schaapcommon::h5parm::JonesParameters;

namespace {
template <size_t PolarizationCount>
void ApplyConjugatedBeam(std::complex<float>* visibilities,
                         const aocommon::MC2x2F& gain1,
                         const aocommon::MC2x2F& gain2);
template <>
void ApplyConjugatedBeam<1>(std::complex<float>* visibilities,
                            const aocommon::MC2x2F& gain1,
                            const aocommon::MC2x2F& gain2) {
  // Stokes-I
  *visibilities = 0.25f * std::conj(aocommon::Trace(gain1)) * (*visibilities) *
                  aocommon::Trace(gain2);
}

template <>
void ApplyConjugatedBeam<4>(std::complex<float>* visibilities,
                            const aocommon::MC2x2F& gain1,
                            const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
  result.AssignTo(visibilities);
}

template <size_t PolarizationCount>
void ApplyBeam(std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
               const aocommon::MC2x2F& gain2);

template <>
void ApplyBeam<1>(std::complex<float>* visibilities,
                  const aocommon::MC2x2F& gain1,
                  const aocommon::MC2x2F& gain2) {
  // Stokes-I
  *visibilities = 0.25f * aocommon::Trace(gain1) * (*visibilities) *
                  std::conj(aocommon::Trace(gain2));
}

template <>
void ApplyBeam<4>(std::complex<float>* visibilities,
                  const aocommon::MC2x2F& gain1,
                  const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.Multiply(visibilities_mc2x2).MultiplyHerm(gain2);
  result.AssignTo(visibilities);
}
}  // namespace

MSGridderBase::~MSGridderBase(){};

MSGridderBase::MSData::MSData()
    : msIndex(0), matchingRows(0), totalRowsProcessed(0) {}

MSGridderBase::MSData::~MSData() {}

MSGridderBase::MSGridderBase(const Settings& settings)
    : _actualInversionWidth(0),
      _actualInversionHeight(0),
      _actualPixelSizeX(0),
      _actualPixelSizeY(0),
      _metaDataCache(nullptr),
      _settings(settings),
      _phaseCentreRA(0.0),
      _phaseCentreDec(0.0),
      _phaseCentreDL(0.0),
      _phaseCentreDM(0.0),
      _facetIndex(0),
      _antennaNames(),
      _imageWidth(0),
      _imageHeight(0),
      _trimWidth(0),
      _trimHeight(0),
      _pixelSizeX(settings.pixelScaleX),
      _pixelSizeY(settings.pixelScaleY),
      _wGridSize(settings.nWLayers),
      _actualWGridSize(0),
      _measurementSets(),
      _dataColumnName(settings.dataColumnName),
      _doImagePSF(false),
      _doSubtractModel(false),
      _addToModel(false),
      _smallInversion(settings.smallInversion),
      _wLimit(settings.wLimit / 100.0),
      _precalculatedWeightInfo(nullptr),
      _polarization(aocommon::Polarization::StokesI),
      _isComplex(false),
      _weighting(settings.weightMode),
      _isFirstIteration(false),
      _visibilityWeightingMode(settings.visibilityWeightingMode),
      _gridMode(GridMode::KaiserBesselKernel),
      _storeImagingWeights(false),
      _theoreticalBeamSize(0.0),
      _hasFrequencies(false),
      _freqHigh(0.0),
      _freqLow(0.0),
      _bandStart(0.0),
      _bandEnd(0.0),
      _startTime(0.0),
      _griddedVisibilityCount(0),
      _totalWeight(0.0),
      _maxGriddedWeight(0.0),
      _visibilityWeightSum(0.0),
      _cachedParmResponse(),
      _h5parm(nullptr),
      _solTabs(std::make_pair(nullptr, nullptr)),
      _h5TimeIndex(std::make_pair(std::numeric_limits<size_t>::max(),
                                  std::numeric_limits<size_t>::max())) {
  computeFacetCentre();
}

void MSGridderBase::setAntennaNames(const casacore::MeasurementSet& ms) {
  _antennaNames.clear();
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    _antennaNames.push_back(antenna.name()(i));
  }
}

int64_t MSGridderBase::getAvailableMemory(double memFraction,
                                          double absMemLimit) {
  static std::atomic<bool> isFirst(true);
  bool printOutput = isFirst.exchange(false);
  long int pageCount = sysconf(_SC_PHYS_PAGES),
           pageSize = sysconf(_SC_PAGE_SIZE);
  int64_t memory = (int64_t)pageCount * (int64_t)pageSize;
  double memSizeInGB = (double)memory / (1024.0 * 1024.0 * 1024.0);
  if (memFraction == 1.0 && absMemLimit == 0.0 && printOutput) {
    Logger::Info << "Detected " << round(memSizeInGB * 10.0) / 10.0
                 << " GB of system memory, usage not limited.\n";
  } else {
    double limitInGB = memSizeInGB * memFraction;
    if (absMemLimit != 0.0 && limitInGB > absMemLimit) limitInGB = absMemLimit;
    if (printOutput) {
      Logger::Info << "Detected " << round(memSizeInGB * 10.0) / 10.0
                   << " GB of system memory, usage limited to "
                   << round(limitInGB * 10.0) / 10.0
                   << " GB (frac=" << round(memFraction * 1000.0) / 10.0
                   << "%, ";
      if (absMemLimit == 0.0)
        Logger::Info << "no limit)\n";
      else
        Logger::Info << "limit=" << round(absMemLimit * 10.0) / 10.0 << "GB)\n";
    }

    memory = int64_t((double)pageCount * (double)pageSize * memFraction);
    if (absMemLimit != 0.0 &&
        double(memory) > double(1024.0 * 1024.0 * 1024.0) * absMemLimit)
      memory = int64_t(double(absMemLimit) * double(1024.0 * 1024.0 * 1024.0));
  }
  return memory;
}

void MSGridderBase::initializeBandData(casacore::MeasurementSet& ms,
                                       MSGridderBase::MSData& msData) {
  msData.bandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
  if (Selection(msData.msIndex).HasChannelRange()) {
    msData.startChannel = Selection(msData.msIndex).ChannelRangeStart();
    msData.endChannel = Selection(msData.msIndex).ChannelRangeEnd();
    Logger::Debug << "Selected channels: " << msData.startChannel << '-'
                  << msData.endChannel << '\n';
    if (msData.startChannel >= msData.bandData.MaxChannels() ||
        msData.endChannel > msData.bandData.MaxChannels() ||
        msData.startChannel == msData.endChannel) {
      std::ostringstream str;
      str << "An invalid channel range was specified! Measurement set only has "
          << msData.bandData.MaxChannels()
          << " channels, requested imaging range is " << msData.startChannel
          << " -- " << msData.endChannel << '.';
      throw std::runtime_error(str.str());
    }
  } else {
    msData.startChannel = 0;
    msData.endChannel = msData.bandData.MaxChannels();
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
  MultiBandData selectedBand = msData.SelectedBand();
  std::vector<float> weightArray(selectedBand.MaxChannels() * NPolInMSProvider);
  double curTimestep = -1, firstTime = -1, lastTime = -1;
  size_t nTimesteps = 0;
  std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);

    if (curTimestep != metaData.time) {
      curTimestep = metaData.time;
      ++nTimesteps;
      if (firstTime == -1) firstTime = curTimestep;
      lastTime = curTimestep;
    }

    const BandData& curBand = selectedBand[metaData.dataDescId];
    double wHi = fabs(metaData.wInM / curBand.SmallestWavelength());
    double wLo = fabs(metaData.wInM / curBand.LongestWavelength());
    double baselineInM =
        sqrt(metaData.uInM * metaData.uInM + metaData.vInM * metaData.vInM +
             metaData.wInM * metaData.wInM);
    double halfWidth = 0.5 * ImageWidth(), halfHeight = 0.5 * ImageHeight();
    if (wHi > msData.maxW || wLo < msData.minW ||
        baselineInM / curBand.SmallestWavelength() > msData.maxBaselineUVW) {
      msReader->ReadWeights(weightArray.data());
      const float* weightPtr = weightArray.data();
      for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
        const double wavelength = curBand.ChannelWavelength(ch);
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
template void MSGridderBase::calculateWLimits<4>(MSGridderBase::MSData& msData);

void MSGridderBase::initializeMSDataVector(
    std::vector<MSGridderBase::MSData>& msDataVector, bool isPredict) {
  if (MeasurementSetCount() == 0)
    throw std::runtime_error(
        "Something is wrong during inversion: no measurement sets given to "
        "inversion algorithm");
  msDataVector = std::vector<MSGridderBase::MSData>(MeasurementSetCount());

  resetMetaData();

  bool hasCache = !_metaDataCache->msDataVector.empty();
  if (!hasCache) _metaDataCache->msDataVector.resize(MeasurementSetCount());
  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    msDataVector[i].msIndex = i;
    initializeMeasurementSet(msDataVector[i], _metaDataCache->msDataVector[i],
                             hasCache, isPredict);
  }

  calculateOverallMetaData(msDataVector.data());
}

void MSGridderBase::initializeMeasurementSet(MSGridderBase::MSData& msData,
                                             MetaDataCache::Entry& cacheEntry,
                                             bool isCacheInitialized,
                                             bool isPredict) {
  MSProvider& msProvider = MeasurementSet(msData.msIndex);
  msData.msProvider = &msProvider;
  SynchronizedMS ms(msProvider.MS());
  if (ms->nrow() == 0) throw std::runtime_error("Table has no rows (no data)");

  initializeBandData(*ms, msData);

  if (HasDenormalPhaseCentre())
    Logger::Debug << "Set has denormal phase centre: dl=" << _phaseCentreDL
                  << ", dm=" << _phaseCentreDM << '\n';

  calculateMSLimits(msData.SelectedBand(), msProvider.StartTime());

  if (isCacheInitialized) {
    msData.maxW = cacheEntry.maxW;
    msData.maxWWithFlags = cacheEntry.maxWWithFlags;
    msData.minW = cacheEntry.minW;
    msData.maxBaselineUVW = cacheEntry.maxBaselineUVW;
    msData.maxBaselineInM = cacheEntry.maxBaselineInM;
    msData.integrationTime = cacheEntry.integrationTime;
  } else {
    if (msProvider.Polarization() == aocommon::Polarization::Instrumental)
      calculateWLimits<4>(msData);
    else
      calculateWLimits<1>(msData);
    cacheEntry.maxW = msData.maxW;
    cacheEntry.maxWWithFlags = msData.maxWWithFlags;
    cacheEntry.minW = msData.minW;
    cacheEntry.maxBaselineUVW = msData.maxBaselineUVW;
    cacheEntry.maxBaselineInM = msData.maxBaselineInM;
    cacheEntry.integrationTime = msData.integrationTime;
  }

#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    // Hard-coded for now
    const bool frequency_interpolation = true;
    const bool use_channel_frequency = true;
    const std::string element_response_string =
        !_settings.beamModel.empty() ? _settings.beamModel : "DEFAULT";

    // Get path to coefficient file in case MWA telescope
    everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
    const std::string coeff_path =
        (telescope_type == everybeam::TelescopeType::kMWATelescope)
            ? wsclean::mwa::FindCoeffFile(_settings.mwaPath)
            : "";

    everybeam::ATermSettings aterm_settings;
    aterm_settings.coeff_path = coeff_path;
    aterm_settings.data_column_name = _settings.dataColumnName;

    everybeam::Options options =
        everybeam::aterms::ATermConfig::ConvertToEBOptions(
            *ms, aterm_settings, frequency_interpolation,
            _settings.useDifferentialLofarBeam, use_channel_frequency,
            element_response_string);

    _telescope = everybeam::Load(*ms, options);
    _pointResponse = _telescope->GetPointResponse(msProvider.StartTime());
    _pointResponse->SetUpdateInterval(_settings.facetBeamUpdateTime);
    _cachedBeamResponse.resize(msData.bandData.MaxChannels() *
                               _pointResponse->GetAllStationsBufferSize());
  } else {
    if (_settings.applyFacetBeam) {
      throw std::runtime_error(
          "-apply-facet-beam was set, but no corresponding facet "
          "regions file was specified.");
    }
    _pointResponse = nullptr;
    _cachedBeamResponse.resize(0);
  }
#else
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    throw std::runtime_error(
        "-apply-facet-beam was set, but wsclean was not compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif

  if (!_settings.facetSolutionFile.empty()) {
    _h5parm.reset(
        new schaapcommon::h5parm::H5Parm(_settings.facetSolutionFile));

    if (_settings.facetSolutionTables.size() == 1) {
      _solTabs = std::make_pair(
          &_h5parm->GetSolTab(_settings.facetSolutionTables[0]), nullptr);
      _correctType =
          JonesParameters::StringToCorrectType(_solTabs.first->GetType());
    } else if (_settings.facetSolutionTables.size() == 2) {
      _solTabs =
          std::make_pair(&_h5parm->GetSolTab(_settings.facetSolutionTables[0]),
                         &_h5parm->GetSolTab(_settings.facetSolutionTables[1]));
      if (_solTabs.first->GetType() != "amplitude") {
        throw std::runtime_error("Type of solution table 0 is " +
                                 _solTabs.first->GetType() +
                                 ", should be 'amplitude'");
      }
      if (_solTabs.second->GetType() != "phase") {
        throw std::runtime_error("Type of solution table 1 is " +
                                 _solTabs.second->GetType() +
                                 ", should be 'phase'");
      }
      _correctType = JonesParameters::CorrectType::FULLJONES;
    } else {
      throw std::runtime_error(
          "Specify the solution table name(s) with "
          "-soltab-names=soltabname1[OPTIONAL,soltabname2]");
    }
  }

  if (isPredict) {
    _degriddingReader = msData.msProvider->MakeReader();
  } else {
    _degriddingReader.reset();
  }
}

void MSGridderBase::calculateOverallMetaData(const MSData* msDataVector) {
  _maxW = 0.0;
  _minW = std::numeric_limits<double>::max();
  double maxBaseline = 0.0;

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    const MSData& msData = msDataVector[i];

    maxBaseline = std::max(maxBaseline, msData.maxBaselineUVW);
    _maxW = std::max(_maxW, msData.maxW);
    _minW = std::min(_minW, msData.minW);
  }
  if (_minW > _maxW) {
    _minW = _maxW;
    Logger::Error
        << "*** Error! ***\n"
           "*** Calculating maximum and minimum w values failed! Make sure the "
           "data selection and scale settings are correct!\n"
           "***\n";
  }

  _theoreticalBeamSize = 1.0 / maxBaseline;
  if (IsFirstIteration()) {
    Logger::Info << "Theoretic beam = "
                 << Angle::ToNiceString(_theoreticalBeamSize) << "\n";
  }
  if (_wLimit != 0.0) {
    _maxW *= (1.0 - _wLimit);
    if (_maxW < _minW) _maxW = _minW;
  }

  if (!HasTrimSize()) SetTrimSize(ImageWidth(), ImageHeight());

  _actualInversionWidth = ImageWidth();
  _actualInversionHeight = ImageHeight();
  _actualPixelSizeX = PixelSizeX();
  _actualPixelSizeY = PixelSizeY();

  if (SmallInversion()) {
    size_t optWidth, optHeight, minWidth, minHeight;
    CalculateFFTSize(_actualInversionWidth, _actualPixelSizeX,
                     _theoreticalBeamSize, minWidth, optWidth);
    CalculateFFTSize(_actualInversionHeight, _actualPixelSizeY,
                     _theoreticalBeamSize, minHeight, optHeight);
    if (optWidth < _actualInversionWidth ||
        optHeight < _actualInversionHeight) {
      size_t newWidth =
          std::max(std::min(optWidth, _actualInversionWidth), size_t(32));
      size_t newHeight =
          std::max(std::min(optHeight, _actualInversionHeight), size_t(32));
      if (IsFirstIteration()) {
        Logger::Info << "Minimal inversion size: " << minWidth << " x "
                     << minHeight << ", using optimal: " << newWidth << " x "
                     << newHeight << "\n";
      }
      _actualPixelSizeX = (double(_actualInversionWidth) * _actualPixelSizeX) /
                          double(newWidth);
      _actualPixelSizeY = (double(_actualInversionHeight) * _actualPixelSizeY) /
                          double(newHeight);
      _actualInversionWidth = newWidth;
      _actualInversionHeight = newHeight;
    } else {
      if (IsFirstIteration()) {
        Logger::Info
            << "Small inversion enabled, but inversion resolution already "
               "smaller than beam size: not using optimization.\n";
      }
    }
  }

  if (IsFirstIteration() || !hasWGridSize()) {
    size_t suggestedGridSize = getSuggestedWGridSize();
    if (!hasWGridSize())
      _actualWGridSize = suggestedGridSize;
    else
      _actualWGridSize = _wGridSize;
  } else
    _actualWGridSize = _wGridSize;
}

template <size_t PolarizationCount>
void MSGridderBase::writeVisibilities(MSProvider& msProvider,
                                      const BandData& curBand,
                                      std::complex<float>* buffer) {
  bool updateRowInH5Parm = true;
#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    updateRowInH5Parm = false;
  }
#endif

  if (_h5parm) {
    if (_antennaNames.empty()) {
      throw std::runtime_error(
          "Antenna names have to be specified in order to apply H5Parm "
          "solutions.");
    }
    MSProvider::MetaData metaData;
    _degriddingReader->ReadMeta(metaData);
    if (updateRowInH5Parm) {
      _degriddingReader->NextInputRow();
    }

    const size_t nparm =
        (_correctType == JonesParameters::CorrectType::FULLJONES) ? 4 : 2;

    // Only update the cached response if one of the time indices in the soltabs
    // changed
    if (_solTabs.first->GetTimeIndex(metaData.time) != _h5TimeIndex.first ||
        (_solTabs.second &&
         _solTabs.second->GetTimeIndex(metaData.time) != _h5TimeIndex.second)) {
      const std::vector<double> freqs(curBand.begin(), curBand.end());
      // FIXME: leads to a small overhead, but this is because _antennaNames are
      // not known in initializeMeasurementSet
      _cachedParmResponse.resize(freqs.size() * _antennaNames.size() * nparm);

      JonesParameters jonesParameters(
          freqs, std::vector<double>{metaData.time}, _antennaNames,
          _correctType, JonesParameters::InterpolationType::NEAREST,
          _facetIndex, _solTabs.first, _solTabs.second, false, 0.0f, 0u,
          JonesParameters::MissingAntennaBehavior::kUnit);
      const auto parms = jonesParameters.GetParms();
      // FIXME: following assignment assumes that the data layout
      // of parms is contiguous in memory, ordered as [station1:pol1:chan1, ...,
      // station1:poln:chan1, ... station2:poln:chan1, stationm:poln:chan1,
      // station2:pol1:chan2, ..., stationn:polm:chank] Check this!
      _cachedParmResponse.assign(&parms(0, 0, 0),
                                 &parms(0, 0, 0) + _cachedParmResponse.size());

      _h5TimeIndex.first = _solTabs.first->GetTimeIndex(metaData.time);
      if (_solTabs.second) {
        _h5TimeIndex.second = _solTabs.second->GetTimeIndex(metaData.time);
      }
    }

    std::complex<float>* iter = buffer;
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      // TODO: template this on nparm for efficiency?
      const size_t offset = ch * _antennaNames.size() * nparm;
      const size_t offset1 = offset + metaData.antenna1 * nparm;
      const size_t offset2 = offset + metaData.antenna2 * nparm;

      if (nparm == 2) {
        ApplyBeam<PolarizationCount>(
            iter,
            aocommon::MC2x2F(_cachedParmResponse[offset1], 0, 0,
                             _cachedParmResponse[offset1 + 1]),
            aocommon::MC2x2F(_cachedParmResponse[offset2], 0, 0,
                             _cachedParmResponse[offset2 + 1]));
      } else {
        ApplyBeam<PolarizationCount>(
            iter, aocommon::MC2x2F(&_cachedParmResponse[offset1]),
            aocommon::MC2x2F(&_cachedParmResponse[offset2]));
      }
      iter += PolarizationCount;
    }
  }

#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    MSProvider::MetaData metaData;
    _degriddingReader->ReadMeta(metaData);
    if (!updateRowInH5Parm) {
      _degriddingReader->NextInputRow();
    }
    _pointResponse->UpdateTime(metaData.time);
    if (_pointResponse->HasTimeUpdate()) {
      if (auto phasedArray =
              dynamic_cast<everybeam::pointresponse::PhasedArrayPoint*>(
                  _pointResponse.get())) {
        phasedArray->UpdateITRFVectors(_facetCentreRA, _facetCentreDec);
      }
      for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
        _pointResponse->CalculateAllStations(
            &_cachedBeamResponse[ch *
                                 _pointResponse->GetAllStationsBufferSize()],
            _facetCentreRA, _facetCentreDec, curBand.ChannelFrequency(ch),
            metaData.fieldId);
      }
    }

    std::complex<float>* iter = buffer;
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
      const size_t offset1 = offset + metaData.antenna1 * 4u;
      const size_t offset2 = offset + metaData.antenna2 * 4u;

      const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
      const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
      ApplyBeam<PolarizationCount>(iter, gain1, gain2);
      iter += PolarizationCount;
    }
  }
#endif

  const bool addToMS = (_facetIndex != 0);
  msProvider.WriteModel(buffer, addToMS);
  msProvider.NextOutputRow();
}

template void MSGridderBase::writeVisibilities<1>(MSProvider& msProvider,
                                                  const BandData& curBand,
                                                  std::complex<float>* buffer);

template void MSGridderBase::writeVisibilities<4>(MSProvider& msProvider,
                                                  const BandData& curBand,
                                                  std::complex<float>* buffer);

template <size_t PolarizationCount>
void MSGridderBase::readAndWeightVisibilities(MSReader& msReader,
                                              InversionRow& rowData,
                                              const BandData& curBand,
                                              float* weightBuffer,
                                              std::complex<float>* modelBuffer,
                                              const bool* isSelected) {
  const std::size_t dataSize = curBand.ChannelCount() * PolarizationCount;
  if (DoImagePSF()) {
    std::fill_n(rowData.data, dataSize, 1.0);
    if (HasDenormalPhaseCentre() && _settings.facetRegionFilename.empty()) {
      const double lmsqrt = std::sqrt(1.0 - PhaseCentreDL() * PhaseCentreDL() -
                                      PhaseCentreDM() * PhaseCentreDM());
      const double shiftFactor = 2.0 * M_PI *
                                 ((rowData.uvw[0] * PhaseCentreDL() +
                                   rowData.uvw[1] * PhaseCentreDM()) +
                                  rowData.uvw[2] * (lmsqrt - 1.0));
      rotateVisibilities<PolarizationCount>(curBand, shiftFactor, rowData.data);
    }
  } else {
    msReader.ReadData(rowData.data);
  }
  rowData.rowId = msReader.RowId();

  if (DoSubtractModel()) {
    msReader.ReadModel(modelBuffer);
    std::complex<float>* modelIter = modelBuffer;
    for (std::complex<float>* iter = rowData.data;
         iter != rowData.data + dataSize; ++iter) {
      *iter -= *modelIter;
      modelIter++;
    }
  }

#ifdef HAVE_EVERYBEAM
  if (_settings.applyFacetBeam && !_settings.facetRegionFilename.empty()) {
    MSProvider::MetaData metaData;
    msReader.ReadMeta(metaData);
    _pointResponse->UpdateTime(metaData.time);
    if (_pointResponse->HasTimeUpdate()) {
      if (auto phasedArray =
              dynamic_cast<everybeam::pointresponse::PhasedArrayPoint*>(
                  _pointResponse.get())) {
        phasedArray->UpdateITRFVectors(_facetCentreRA, _facetCentreDec);
      }
      for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
        _pointResponse->CalculateAllStations(
            &_cachedBeamResponse[ch *
                                 _pointResponse->GetAllStationsBufferSize()],
            _facetCentreRA, _facetCentreDec, curBand.ChannelFrequency(ch),
            metaData.fieldId);
      }
    }

    // rowData.data contains the visibilities
    std::complex<float>* iter = rowData.data;
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
      const size_t offset1 = offset + metaData.antenna1 * 4u;
      const size_t offset2 = offset + metaData.antenna2 * 4u;

      const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
      const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
      ApplyConjugatedBeam<PolarizationCount>(iter, gain1, gain2);
      iter += PolarizationCount;
    }
  }
#endif

  if (_h5parm) {
    if (_antennaNames.empty()) {
      throw std::runtime_error(
          "Antenna names have to be specified in order to apply H5Parm "
          "solutions.");
    }
    MSProvider::MetaData metaData;
    msReader.ReadMeta(metaData);

    const size_t nparm =
        (_correctType == JonesParameters::CorrectType::FULLJONES) ? 4 : 2;

    // Only update the cached response if one of the time indices in the soltabs
    // changed
    if (_solTabs.first->GetTimeIndex(metaData.time) != _h5TimeIndex.first ||
        (_solTabs.second &&
         _solTabs.second->GetTimeIndex(metaData.time) != _h5TimeIndex.second)) {
      const std::vector<double> freqs(curBand.begin(), curBand.end());
      // FIXME: leads to a small overhead, but this is because _antennaNames are
      // not known in initializeMeasurementSet
      _cachedParmResponse.resize(freqs.size() * _antennaNames.size() * nparm);

      JonesParameters jonesParameters(
          freqs, std::vector<double>{metaData.time}, _antennaNames,
          _correctType, JonesParameters::InterpolationType::NEAREST,
          _facetIndex, _solTabs.first, _solTabs.second, false, 0.0f, 0u,
          JonesParameters::MissingAntennaBehavior::kUnit);
      const auto parms = jonesParameters.GetParms();
      // FIXME: following assignment assumes that the data layout
      // of parms is contiguous in mem, ordered as [station1:pol1:chan1, ...,
      // station1:poln:chan1, ... station2:poln:chan1, stationm:poln:chan1,
      // station2:pol1:chan2, ..., stationn:polm:chank] Check this!
      _cachedParmResponse.assign(&parms(0, 0, 0),
                                 &parms(0, 0, 0) + _cachedParmResponse.size());

      _h5TimeIndex.first = _solTabs.first->GetTimeIndex(metaData.time);
      if (_solTabs.second) {
        _h5TimeIndex.second = _solTabs.second->GetTimeIndex(metaData.time);
      }
    }

    std::complex<float>* iter = rowData.data;
    for (size_t ch = 0; ch < curBand.ChannelCount(); ++ch) {
      // TODO: template this on nparm for efficiency?
      aocommon::MC2x2F gain1;
      aocommon::MC2x2F gain2;
      const size_t offset = ch * _antennaNames.size() * nparm;
      const size_t offset1 = offset + metaData.antenna1 * nparm;
      const size_t offset2 = offset + metaData.antenna2 * nparm;

      if (nparm == 2) {
        ApplyConjugatedBeam<PolarizationCount>(
            iter,
            aocommon::MC2x2F(_cachedParmResponse[offset1], 0, 0,
                             _cachedParmResponse[offset1 + 1]),
            aocommon::MC2x2F(_cachedParmResponse[offset2], 0, 0,
                             _cachedParmResponse[offset2 + 1]));
      } else {
        ApplyConjugatedBeam<PolarizationCount>(
            iter, aocommon::MC2x2F(&_cachedParmResponse[offset1]),
            aocommon::MC2x2F(&_cachedParmResponse[offset2]));
      }
      iter += PolarizationCount;
    }
  }

  msReader.ReadWeights(weightBuffer);

  // Any visibilities that are not gridded in this pass
  // should not contribute to the weight sum, so set these
  // to have zero weight.
  for (size_t ch = 0; ch != dataSize; ++ch) {
    if (!isSelected[ch]) weightBuffer[ch] = 0.0;
  }

  switch (VisibilityWeightingMode()) {
    case VisibilityWeightingMode::NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case VisibilityWeightingMode::SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t chp = 0; chp != dataSize; ++chp)
        weightBuffer[chp] *= weightBuffer[chp];
      break;
    case VisibilityWeightingMode::UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t chp = 0; chp != dataSize; ++chp) {
        if (weightBuffer[chp] != 0.0) weightBuffer[chp] = 1.0f;
      }
      break;
  }

  // Calculate imaging weights
  std::complex<float>* dataIter = rowData.data;
  float* weightIter = weightBuffer;
  _scratchWeights.resize(curBand.ChannelCount());
  for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
    double u = rowData.uvw[0] / curBand.ChannelWavelength(ch),
           v = rowData.uvw[1] / curBand.ChannelWavelength(ch),
           imageWeight = GetImageWeights()->GetWeight(u, v);
    _scratchWeights[ch] = imageWeight;

    for (size_t p = 0; p != PolarizationCount; ++p) {
      double cumWeight = *weightIter * imageWeight;
      if (p == 0 && cumWeight != 0.0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        _visibilityWeightSum += *weightIter;
        _maxGriddedWeight = std::max(cumWeight, _maxGriddedWeight);
        ++_griddedVisibilityCount;
        // Total weight includes imaging weights
        _totalWeight += cumWeight;
      }
      *weightIter = cumWeight;
      *dataIter *= *weightIter;
      ++dataIter;
      ++weightIter;
    }
  }
  if (StoreImagingWeights())
    msReader.WriteImagingWeights(_scratchWeights.data());
}

template void MSGridderBase::readAndWeightVisibilities<1>(
    MSReader& msReader, InversionRow& newItem, const BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<4>(
    MSReader& msReader, InversionRow& newItem, const BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected);

template <size_t PolarizationCount>
void MSGridderBase::rotateVisibilities(const BandData& bandData,
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
    const BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<4>(
    const BandData& bandData, double shiftFactor,
    std::complex<float>* dataIter);