#include "msgridderbase.h"

#include "msgriddermanager.h"
#include "msprovidercollection.h"

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

namespace wsclean {

// Defined out of class to allow the class the be used in a std::unique_ptr.
MSGridderBase::~MSGridderBase() = default;

MSGridderBase::MSGridderBase(const Settings& settings,
                             MsProviderCollection& ms_provider_collection)
    : ms_data_vector_(ms_provider_collection.ms_data_vector_),
      settings_(settings),
      w_grid_size_(settings.nWLayers),
      data_column_name_(settings.dataColumnName),
      small_inversion_(settings.minGridResolution),
      weighting_(settings.weightMode),
      visibility_weighting_mode_(settings.visibilityWeightingMode) {
#ifdef HAVE_EVERYBEAM
  visibility_modifier_.SetBeamInfo(settings.beamMode,
                                   settings.beamNormalisationMode);
#endif
}

void MSGridderBase::initializePointResponse(
    const MsProviderCollection::MsData& msData) {
#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    const std::string element_response_string =
        !settings_.beamModel.empty() ? settings_.beamModel : "DEFAULT";
    visibility_modifier_.InitializePointResponse(
        msData.ms_provider->MS(), settings_.facetBeamUpdateTime,
        element_response_string, msData.band_data.ChannelCount(),
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

void MSGridderBase::StartMeasurementSet(
    const MsProviderCollection::MsData& msData, bool isPredict) {
  initializePointResponse(msData);
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.InitializeCacheParmResponse(
        msData.antenna_names, msData.SelectedBand(), msData.original_ms_index);
  }
  original_ms_index_ = msData.original_ms_index;
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

void MSGridderBase::ResetVisibilityModifierCache() {
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.ResetCache(GetMsCount());
  }
  visibility_modifier_.ResetSums();
}

void MSGridderBase::calculateOverallMetaData() {
  theoretical_beam_size_ = 1.0 / MaxBaseline();
  if (IsFirstTask()) {
    Logger::Info << "Theoretic beam = " +
                        aocommon::units::Angle::ToNiceString(
                            theoretical_beam_size_) +
                        "\n";
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
        Logger::Info << "Minimal inversion size: " + std::to_string(minWidth) +
                            " x " + std::to_string(minHeight) +
                            ", using optimal: " + std::to_string(newWidth) +
                            " x " + std::to_string(newHeight) + "\n";
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

template <GainMode Mode>
void MSGridderBase::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData) {
  assert(GetPsfMode() == PsfMode::kNone);  // The PSF is never predicted.

#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam) {
    visibility_modifier_.CacheBeamResponse(metaData.time, metaData.fieldId,
                                           curBand);

    visibility_modifier_.ApplyBeamResponse<Mode>(
        buffer, curBand.ChannelCount(), metaData.antenna1, metaData.antenna2);
  }
#endif

  if (visibility_modifier_.HasH5Parm()) {
    assert(!settings_.facetRegionFilename.empty());
    visibility_modifier_.CacheParmResponse(metaData.time, curBand,
                                           original_ms_index_);

    visibility_modifier_.ApplyParmResponse<Mode>(
        buffer, original_ms_index_, curBand.ChannelCount(), n_antennas,
        metaData.antenna1, metaData.antenna2);
  }

  {
    const size_t lock_index =
        facet_group_index_ * GetMsCount() + original_ms_index_;
    std::unique_ptr<GriddingTaskManager::WriterLock> lock =
        writer_lock_manager_->GetLock(lock_index);
    ms_provider.WriteModel(buffer, IsFacet());
  }
  ms_provider.NextOutputRow();
}

template void MSGridderBase::WriteInstrumentalVisibilities<GainMode::kXX>(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void MSGridderBase::WriteInstrumentalVisibilities<GainMode::kYY>(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void MSGridderBase::WriteInstrumentalVisibilities<GainMode::kTrace>(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void
MSGridderBase::WriteInstrumentalVisibilities<GainMode::k2VisDiagonal>(
    MSProvider& ms_provider, size_t n_antennas,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    MSProvider::MetaData& metaData);

template void MSGridderBase::WriteInstrumentalVisibilities<GainMode::kFull>(
    MSProvider& ms_provider, size_t n_antennas,
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
      // The point source is shifted to the centre of the current DdPsf
      // position
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

}  // namespace wsclean
