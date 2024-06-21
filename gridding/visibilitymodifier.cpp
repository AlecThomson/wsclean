#include "visibilitymodifier.h"

#include "../msproviders/synchronizedms.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/scalar/vector4.h>

using aocommon::HMatrix4x4;
using aocommon::MC2x2;
using aocommon::MC2x2F;

namespace {
void setNonFiniteToZero(std::vector<std::complex<float>>& values) {
  for (std::complex<float>& v : values) {
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
      v = 0.0;
    }
  }
}

constexpr bool ShouldApplyCorrection(ModifierBehaviour behaviour) {
  return behaviour == ModifierBehaviour::kApply ||
         behaviour == ModifierBehaviour::kApplyAndSum;
}

constexpr bool ShouldSumCorrection(ModifierBehaviour behaviour) {
  return behaviour == ModifierBehaviour::kSum ||
         behaviour == ModifierBehaviour::kApplyAndSum;
}

/**
 * @brief Apply gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 2 or 4 for IDG, 1 for all other
 * gridders.
 * @tparam Mode Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <GainMode Mode>
void ApplyGain(std::complex<float>* visibilities, const MC2x2F& gain1,
               const MC2x2F& gain2) {
  if constexpr (Mode == GainMode::kXX) {
    *visibilities = gain1[0] * (*visibilities) * std::conj(gain2[0]);
  } else if constexpr (Mode == GainMode::kYY) {
    *visibilities = gain1[3] * (*visibilities) * std::conj(gain2[3]);
  } else if constexpr (Mode == GainMode::kTrace) {
    // Stokes-I. Have to calculate v' = G1 x v x G2^H:
    // v' = 0.5 (V_xx + V_yy) with V = v x (G1 x G2^H)
    // V_xx = v x (g1_xx g2_xx* + g1_yx g2_yx*), V_yy = v x (g1_xy g2_xy* +
    // g1_yy g2_yy*). Hence v' = 0.5 * double_dot(G1, G2*)
    *visibilities *= 0.5f * gain1.DoubleDot(gain2.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = gain1[0] * visibilities[0] * std::conj(gain2[0]) +
                      gain1[1] * visibilities[1] * std::conj(gain2[1]);
    visibilities[1] = gain1[3] * visibilities[1] * std::conj(gain2[3]) +
                      gain1[2] * visibilities[0] * std::conj(gain2[2]);
  } else if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    const MC2x2F visibilities_mc2x2(visibilities);
    const MC2x2F result =
        gain1.Multiply(visibilities_mc2x2).MultiplyHerm(gain2);
    result.AssignTo(visibilities);
  }
}

/**
 * @brief Apply conjugated gains to the visibilities.
 *
 * @tparam Mode Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <GainMode Mode>
void ApplyConjugatedGain(std::complex<float>* visibilities, const MC2x2F& gain1,
                         const MC2x2F& gain2) {
  if constexpr (Mode == GainMode::kXX) {
    *visibilities = std::conj(gain1[0]) * (*visibilities) * gain2[0];
  } else if constexpr (Mode == GainMode::kYY) {
    *visibilities = std::conj(gain1[3]) * (*visibilities) * gain2[3];
  } else if constexpr (Mode == GainMode::kTrace) {
    // See calculation in ApplyGain() for explanation of double dot.
    *visibilities *= 0.5f * gain2.DoubleDot(gain1.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = std::conj(gain1[0]) * visibilities[0] * gain2[0] +
                      std::conj(gain1[2]) * visibilities[1] * gain2[2];
    visibilities[1] = std::conj(gain1[3]) * visibilities[1] * gain2[3] +
                      std::conj(gain1[1]) * visibilities[0] * gain2[1];
  } else if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    const MC2x2F visibilities_mc2x2(visibilities);
    const MC2x2F result =
        gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
    result.AssignTo(visibilities);
  }
}

}  // namespace

void VisibilityModifier::InitializePointResponse(
    SynchronizedMS&& ms, double facet_beam_update_time,
    const std::string& element_response, size_t n_channels,
    const std::string& data_column, const std::string& mwa_path) {
#ifdef HAVE_EVERYBEAM
  // Hard-coded for now
  const bool frequency_interpolation = true;
  const bool use_channel_frequency = true;

  // Get path to coefficient file for MWA telescope
  everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
  const std::string coeff_path =
      (telescope_type == everybeam::TelescopeType::kMWATelescope)
          ? wsclean::mwa::FindCoeffFile(mwa_path)
          : "";

  everybeam::ATermSettings aterm_settings;
  aterm_settings.coeff_path = coeff_path;
  aterm_settings.data_column_name = data_column;

  everybeam::Options options =
      everybeam::aterms::ATermConfig::ConvertToEBOptions(
          *ms, aterm_settings, frequency_interpolation, _beamNormalisationMode,
          use_channel_frequency, element_response, _beamModeString);
  _beamMode = options.beam_mode;

  _telescope = everybeam::Load(*ms, options);
  // Initialize with 0.0 time to make sure first call to UpdateTime()
  // will fill the beam response cache.
  _pointResponse = _telescope->GetPointResponse(0.0);
  _pointResponse->SetUpdateInterval(facet_beam_update_time);
  _pointResponseBufferSize = _pointResponse->GetAllStationsBufferSize();
  _cachedBeamResponse.resize(n_channels * _pointResponseBufferSize);
#endif
}

void VisibilityModifier::InitializeMockResponse(
    size_t n_channels, size_t n_stations,
    const std::vector<std::complex<double>>& beam_response,
    const std::vector<std::complex<float>>& parm_response) {
  _pointResponseBufferSize = n_stations * 4;
  assert(beam_response.size() == n_channels * _pointResponseBufferSize);
  assert(parm_response.size() == n_channels * n_stations * 2 ||
         parm_response.size() == n_channels * n_stations * 4);
#ifdef HAVE_EVERYBEAM
  _cachedBeamResponse.assign(beam_response.begin(), beam_response.end());
#endif
  _cachedParmResponse.emplace(0, parm_response);
  _timeOffsets = {std::pair(0, 0)};
}

void VisibilityModifier::CacheParmResponse(
    double time, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& band, size_t ms_index) {
  using schaapcommon::h5parm::JonesParameters;

  const size_t solution_index = (*_h5parms).size() == 1 ? 0 : ms_index;

  // Only extract DD solutions if the corresponding cache entry is empty.
  if (_cachedParmResponse[ms_index].empty()) {
    const size_t nparms = NValuesPerSolution(ms_index);
    const std::vector<double> freqs(band.begin(), band.end());
    const size_t responseSize = _cachedMSTimes[ms_index].size() * freqs.size() *
                                antennaNames.size() * nparms;
    const std::string dirName = (*_h5parms)[solution_index].GetNearestSource(
        _facetDirectionRA, _facetDirectionDec);
    schaapcommon::h5parm::SolTab* const first_solution =
        (*_firstSolutions)[solution_index];
    schaapcommon::h5parm::SolTab* const second_solution =
        _secondSolutions->empty() ? nullptr
                                  : (*_secondSolutions)[solution_index];
    const size_t dirIndex = first_solution->GetDirIndex(dirName);
    JonesParameters jonesParameters(
        freqs, _cachedMSTimes[ms_index], antennaNames,
        (*_gainTypes)[solution_index],
        JonesParameters::InterpolationType::NEAREST, dirIndex, first_solution,
        second_solution, false, 0.0f, 0u,
        JonesParameters::MissingAntennaBehavior::kUnit);
    // parms (Casacore::Cube) is column major
    const casacore::Cube<std::complex<float>>& parms =
        jonesParameters.GetParms();
    _cachedParmResponse[ms_index].assign(&parms(0, 0, 0),
                                         &parms(0, 0, 0) + responseSize);
    setNonFiniteToZero(_cachedParmResponse[ms_index]);
  }

  const auto it =
      std::find(_cachedMSTimes[ms_index].begin() + _timeOffsets[ms_index],
                _cachedMSTimes[ms_index].end(), time);
  if (it != _cachedMSTimes[ms_index].end()) {
    // Update _timeOffset value with index
    _timeOffsets[ms_index] =
        std::distance(_cachedMSTimes[ms_index].begin(), it);
  } else {
    throw std::runtime_error(
        "Time not found in cached times. A potential reason could be that the "
        "time values in the provided MS are not in ascending order. "
        "Error occurred with ms index = " +
        std::to_string(ms_index) +
        ", cache "
        "contained " +
        std::to_string(_cachedMSTimes[ms_index].size()) + " elements.\n");
  }
}

#ifdef HAVE_EVERYBEAM
void VisibilityModifier::CacheBeamResponse(double time, size_t fieldId,
                                           const aocommon::BandData& band) {
  _pointResponse->UpdateTime(time);
  if (_pointResponse->HasTimeUpdate()) {
    for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
      _pointResponse->ResponseAllStations(
          _beamMode, &_cachedBeamResponse[ch * _pointResponseBufferSize],
          _facetDirectionRA, _facetDirectionDec, band.ChannelFrequency(ch),
          fieldId);
    }
  }
}

template <GainMode Mode>
void VisibilityModifier::ApplyBeamResponse(std::complex<float>* data,
                                           size_t n_channels, size_t antenna1,
                                           size_t antenna2) {
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * _pointResponseBufferSize;
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const MC2x2F gain2(&_cachedBeamResponse[offset2]);
    ApplyGain<Mode>(data, gain1, gain2);
    data += GetNVisibilities(Mode);
  }
}

template void VisibilityModifier::ApplyBeamResponse<GainMode::kXX>(
    std::complex<float>* data, size_t n_channels, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kYY>(
    std::complex<float>* data, size_t n_channels, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kTrace>(
    std::complex<float>* data, size_t n_channels, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, size_t n_channels, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kFull>(
    std::complex<float>* data, size_t n_channels, size_t antenna1,
    size_t antenna2);

template <ModifierBehaviour Behaviour, GainMode Mode>
void VisibilityModifier::ApplyConjugatedBeamResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward) {
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * _pointResponseBufferSize;
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const MC2x2F gain2(&_cachedBeamResponse[offset2]);
    if constexpr (ShouldApplyCorrection(Behaviour)) {
      if (apply_forward) {
        ApplyGain<Mode>(data, gain1, gain2);
      }
      ApplyConjugatedGain<Mode>(data, gain1, gain2);
    }
    if constexpr (ShouldSumCorrection(Behaviour)) {
      // This assumes that the weights of the polarizations are the same
      correction_sum_.Add<Mode>(gain1, gain2, image_weights[ch] * weights[0]);
    }

    data += GetNVisibilities(Mode);
    weights += GetNVisibilities(Mode);
  }
}

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApply, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApply, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApply, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kSum, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApply, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kSum, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApply, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kSum, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);
#endif

template <GainMode Mode>
void VisibilityModifier::ApplyParmResponse(std::complex<float>* data,
                                           size_t ms_index, size_t n_channels,
                                           size_t n_antennas, size_t antenna1,
                                           size_t antenna2) {
  const size_t nparms = NValuesPerSolution(ms_index);
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffsets[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(_cachedParmResponse[ms_index][offset1], 0, 0,
                         _cachedParmResponse[ms_index][offset1 + 1]);
      const MC2x2F gain2(_cachedParmResponse[ms_index][offset2], 0, 0,
                         _cachedParmResponse[ms_index][offset2 + 1]);
      ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffsets[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(&_cachedParmResponse[ms_index][offset1]);
      const MC2x2F gain2(&_cachedParmResponse[ms_index][offset2]);
      ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
    }
  }
}

template void VisibilityModifier::ApplyParmResponse<GainMode::kXX>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyParmResponse<GainMode::kYY>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyParmResponse<GainMode::kTrace>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyParmResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyParmResponse<GainMode::kFull>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template <ModifierBehaviour Behaviour, GainMode Mode>
void VisibilityModifier::ApplyConjugatedParmResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward) {
  const size_t nparms = NValuesPerSolution(ms_index);

  // Conditional could be templated once C++ supports partial function
  // specialization
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffsets[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(_cachedParmResponse[ms_index][offset1], 0, 0,
                         _cachedParmResponse[ms_index][offset1 + 1]);
      const MC2x2F gain2(_cachedParmResponse[ms_index][offset2], 0, 0,
                         _cachedParmResponse[ms_index][offset2 + 1]);
      if constexpr (ShouldApplyCorrection(Behaviour)) {
        if (apply_forward) {
          ApplyGain<Mode>(data, gain1, gain2);
        }
        ApplyConjugatedGain<Mode>(data, gain1, gain2);
      }
      if constexpr (ShouldSumCorrection(Behaviour)) {
        // This multiplies a lot of zeros so could be done more efficiently
        correction_sum_.Add<Mode>(gain1, gain2, image_weights[ch] * weights[0]);
      }

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffsets[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(&_cachedParmResponse[ms_index][offset1]);
      const MC2x2F gain2(&_cachedParmResponse[ms_index][offset2]);
      if constexpr (ShouldApplyCorrection(Behaviour)) {
        if (apply_forward) {
          ApplyGain<Mode>(data, gain1, gain2);
        }
        ApplyConjugatedGain<Mode>(data, gain1, gain2);
      }
      if constexpr (ShouldSumCorrection(Behaviour)) {
        // Assumes that the weights of the polarizations are the same
        correction_sum_.Add<Mode>(gain1, gain2, image_weights[ch] * weights[0]);
      }

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  }
}

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApply, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApply, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApply, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kSum, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApply, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kSum, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApply, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kSum, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<
    ModifierBehaviour::kApplyAndSum, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

#ifdef HAVE_EVERYBEAM
template <ModifierBehaviour Behaviour, GainMode Mode>
void VisibilityModifier::ApplyConjugatedDual(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward) {
  const size_t nparms = NValuesPerSolution(ms_index);

  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Compute facet beam
      const size_t beam_offset = ch * _pointResponseBufferSize;
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);

      // Get H5 solutions
      // Column major indexing
      const size_t h5_offset =
          (_timeOffsets[ms_index] * n_channels + ch) * n_stations * nparms;
      const size_t h5_offset1 = h5_offset + antenna1 * nparms;
      const size_t h5_offset2 = h5_offset + antenna2 * nparms;
      const MC2x2F gain_h5_1(_cachedParmResponse[ms_index][h5_offset1], 0, 0,
                             _cachedParmResponse[ms_index][h5_offset1 + 1]);
      const MC2x2F gain_h5_2(_cachedParmResponse[ms_index][h5_offset2], 0, 0,
                             _cachedParmResponse[ms_index][h5_offset2 + 1]);

      // Combine H5parm and beam. The beam is applied first on the data,
      // and therefore needs to be the last in the multiplication.
      const MC2x2F gain_combined_1 = gain_h5_1 * gain_b_1;
      const MC2x2F gain_combined_2 = gain_h5_2 * gain_b_2;

      if constexpr (ShouldApplyCorrection(Behaviour)) {
        if (apply_forward) {
          ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
        }
        ApplyConjugatedGain<Mode>(data, gain_combined_1, gain_combined_2);
      }
      if constexpr (ShouldSumCorrection(Behaviour)) {
        beam_correction_sum_.Add<Mode>(gain_b_1, gain_b_2,
                                       weights[0] * image_weights[ch]);
        correction_sum_.Add<Mode>(gain_combined_1, gain_combined_2,
                                  weights[0] * image_weights[ch]);
      }

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  } else {
    // This branch handles full jones H5 parm files (nparms == 4)
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Get facet beam
      const size_t beam_offset = ch * _pointResponseBufferSize;
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);

      // Get h5 solution
      // Column major indexing
      const size_t offset_h5 =
          (_timeOffsets[ms_index] * n_channels + ch) * n_stations * nparms;
      const size_t offset_h5_1 = offset_h5 + antenna1 * nparms;
      const size_t offset_h5_2 = offset_h5 + antenna2 * nparms;
      const MC2x2F gain_h5_1(&_cachedParmResponse[ms_index][offset_h5_1]);
      const MC2x2F gain_h5_2(&_cachedParmResponse[ms_index][offset_h5_2]);

      // Combine H5parm and beam. The beam is applied first on the data,
      // and therefore needs to be the last in the multiplication.
      const MC2x2F gain_combined_1 = gain_h5_1 * gain_b_1;
      const MC2x2F gain_combined_2 = gain_h5_2 * gain_b_2;

      if constexpr (ShouldApplyCorrection(Behaviour)) {
        ApplyConjugatedGain<Mode>(data, gain_combined_1, gain_combined_2);
        if (apply_forward) {
          ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
        }
      }
      if constexpr (ShouldSumCorrection(Behaviour)) {
        beam_correction_sum_.Add<Mode>(gain_b_1, gain_b_2,
                                       weights[0] * image_weights[ch]);
        correction_sum_.Add<Mode>(gain_combined_1, gain_combined_2,
                                  weights[0] * image_weights[ch]);
      }
      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  }
}

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kApply,
                                                      GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void
VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<
    ModifierBehaviour::kApplyAndSum, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kApply,
                                                      GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void
VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<
    ModifierBehaviour::kApplyAndSum, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kApply,
                                                      GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kSum,
                                                      GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<
    ModifierBehaviour::kApplyAndSum, GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kApply,
                                                      GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kSum,
                                                      GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<
    ModifierBehaviour::kApplyAndSum, GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kApply,
                                                      GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<ModifierBehaviour::kSum,
                                                      GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<
    ModifierBehaviour::kApplyAndSum, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

#endif  // HAVE_EVERYBEAM
