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

using aocommon::MC2x2F;

namespace {
void setNonFiniteToZero(std::vector<std::complex<float>>& values) {
  for (std::complex<float>& v : values) {
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
      v = 0.0;
    }
  }
}

/**
 * @brief Compute the gain from the given solution matrices.
 *
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account? See GainMode for further documentation.
 */
template <GainMode Mode>
std::complex<float> ComputeGain(const MC2x2F& gain1, const MC2x2F& gain2) {
  if constexpr (Mode == GainMode::kXX)
    return gain2[0] * std::conj(gain1[0]);
  else if constexpr (Mode == GainMode::kYY)
    return gain2[3] * std::conj(gain1[3]);
  else if constexpr (Mode == GainMode::kTrace ||
                     Mode == GainMode::k2VisDiagonal)
    return 0.5f *
           (gain2[0] * std::conj(gain1[0]) + gain2[3] * std::conj(gain1[3]));
  else
    throw std::runtime_error("ComputeGain<GainMode::kFull> not implemented!");
}

/**
 * @brief Compute the weighted squared gain based on the given gain matrices.
 *
 * @tparam Mode Which entry or entries from the gain matrices should be
 * taken into account? See GainMode for further documentation.
 */
template <GainMode Mode>
float ComputeWeightedSquaredGain(const MC2x2F& gain1, const MC2x2F& gain2,
                                 const float* weights) {
  if constexpr (Mode == GainMode::k2VisDiagonal) {
    // The real gain is the average of these two terms, so requires an
    // extra factor of 2. This is however corrected in the normalization
    // later on.
    return std::norm(gain1[0]) * std::norm(gain2[0]) * weights[0] +
           std::norm(gain1[3]) * std::norm(gain2[3]) * weights[1];
    // TODO is this better? :
    // const std::complex<float> g = ComputeGain<Mode>(gain1, gain2);
    // return std::norm(g) * (weights[0] + weights[1]);
  } else {  // XX / YY / trace
    const std::complex<float> g = ComputeGain<Mode>(gain1, gain2);
    return std::norm(g) * *weights;
  }
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
    // Stokes-I. Have to calculate v <- G1 x v x G2^H:
    // v <- 0.5 (V_xx + V_yy) with V = v x (G1 x G2^H)
    // V_xx = v x (g1_xx g2_xx* + g1_yx g2_yx*), V_yy = v x (g1_xy g2_xy* +
    // g1_yy g2_yy*) Hence v <- 0.5 * double_dot(G1, G2*)
    *visibilities *= 0.5f * gain1.DoubleDot(gain2.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = gain1[0] * visibilities[0] * std::conj(gain2[0]);
    visibilities[1] = gain1[3] * visibilities[1] * std::conj(gain2[3]);
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
 * @tparam PolarizationCount polarization count, 2 or 4 for IDG, 1 for all other
 * gridders.
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
    *visibilities *= 0.5f * gain2.DoubleDot(gain1.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = std::conj(gain1[0]) * visibilities[0] * gain2[0];
    visibilities[1] = std::conj(gain1[3]) * visibilities[1] * gain2[3];
  } else if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    const MC2x2F visibilities_mc2x2(visibilities);
    const MC2x2F result =
        gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
    result.AssignTo(visibilities);
  }
}

/**
 * @brief Apply conjugated gains to the visibilities.
 *
 * @tparam GainEntry decides which entry or entries from the gain matrices
 * should be taken into account in the product, only the diagonal, or the full
 * matrix? See also the documentation of GainMode.
 */
template <GainMode Mode>
MC2x2F MultiplyGains(const MC2x2F& gain_a, const MC2x2F& gain_b) {
  if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    return gain_a.Multiply(gain_b);
  } else {
    return MC2x2F(gain_a[0] * gain_b[0], 0, 0, gain_a[3] * gain_b[3]);
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
  _cachedParmResponse.emplace_back(parm_response);
  _timeOffset = {0};
}

void VisibilityModifier::CacheParmResponse(
    double time, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& band, size_t ms_index) {
  using schaapcommon::h5parm::JonesParameters;

  const size_t solution_index = (*_h5parms).size() == 1 ? 0 : ms_index;

  // Only extract DD solutions if the corresponding cache entry is empty.
  if (_cachedParmResponse[ms_index].empty()) {
    const size_t nparms = NParms(ms_index);
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
      std::find(_cachedMSTimes[ms_index].begin() + _timeOffset[ms_index],
                _cachedMSTimes[ms_index].end(), time);
  if (it != _cachedMSTimes[ms_index].end()) {
    // Update _timeOffset value with index
    _timeOffset[ms_index] = std::distance(_cachedMSTimes[ms_index].begin(), it);
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

template <GainMode Mode>
void VisibilityModifier::ApplyConjugatedBeamResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward) {
  double local_correction_sum = 0.0;
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * _pointResponseBufferSize;
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const MC2x2F gain2(&_cachedBeamResponse[offset2]);
    if (apply_forward) {
      ApplyGain<Mode>(data, gain1, gain2);
    }
    ApplyConjugatedGain<Mode>(data, gain1, gain2);
    const float weighted_squared_gain =
        ComputeWeightedSquaredGain<Mode>(gain1, gain2, weights);
    local_correction_sum += image_weights[ch] * weighted_squared_gain;

    data += GetNVisibilities(Mode);
    weights += GetNVisibilities(Mode);
  }
  correction_sum_ += local_correction_sum;
}

template void VisibilityModifier::ApplyConjugatedBeamResponse<GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void
VisibilityModifier::ApplyConjugatedBeamResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedBeamResponse<GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward);

#endif

template <GainMode Mode>
void VisibilityModifier::ApplyParmResponse(std::complex<float>* data,
                                           size_t ms_index, size_t n_channels,
                                           size_t n_antennas, size_t antenna1,
                                           size_t antenna2) {
  const size_t nparms = NParms(ms_index);
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * n_channels + ch) * n_antennas * nparms;
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
          (_timeOffset[ms_index] * n_channels + ch) * n_antennas * nparms;
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

template <GainMode Mode>
void VisibilityModifier::ApplyConjugatedParmResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward) {
  const size_t nparms = NParms(ms_index);

  double local_correction_sum = 0.0;

  // Conditional could be templated once C++ supports partial function
  // specialization
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(_cachedParmResponse[ms_index][offset1], 0, 0,
                         _cachedParmResponse[ms_index][offset1 + 1]);
      const MC2x2F gain2(_cachedParmResponse[ms_index][offset2], 0, 0,
                         _cachedParmResponse[ms_index][offset2 + 1]);
      if (apply_forward) {
        ApplyGain<Mode>(data, gain1, gain2);
      }
      ApplyConjugatedGain<Mode>(data, gain1, gain2);
      const float weighted_squared_gain =
          ComputeWeightedSquaredGain<Mode>(gain1, gain2, weights);
      local_correction_sum += image_weights[ch] * weighted_squared_gain;

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(&_cachedParmResponse[ms_index][offset1]);
      const MC2x2F gain2(&_cachedParmResponse[ms_index][offset2]);
      if (apply_forward) {
        ApplyGain<Mode>(data, gain1, gain2);
      }
      ApplyConjugatedGain<Mode>(data, gain1, gain2);
      const float weighted_squared_gain =
          ComputeWeightedSquaredGain<Mode>(gain1, gain2, weights);
      local_correction_sum += image_weights[ch] * weighted_squared_gain;

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  }
  correction_sum_ += local_correction_sum;
}

template void VisibilityModifier::ApplyConjugatedParmResponse<GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void
VisibilityModifier::ApplyConjugatedParmResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedParmResponse<GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward);

#ifdef HAVE_EVERYBEAM
template <GainMode Mode>
void VisibilityModifier::ApplyConjugatedDual(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward) {
  double correctionSum = 0.0;
  double h5Sum = 0.0;
  const size_t nparms = NParms(ms_index);

  // Conditional could be templated once C++ supports partial function
  // specialization
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
          (_timeOffset[ms_index] * n_channels + ch) * n_stations * nparms;
      const size_t h5_offset1 = h5_offset + antenna1 * nparms;
      const size_t h5_offset2 = h5_offset + antenna2 * nparms;
      const MC2x2F gain_h5_1(_cachedParmResponse[ms_index][h5_offset1], 0, 0,
                             _cachedParmResponse[ms_index][h5_offset1 + 1]);
      const MC2x2F gain_h5_2(_cachedParmResponse[ms_index][h5_offset2], 0, 0,
                             _cachedParmResponse[ms_index][h5_offset2 + 1]);

      // Combine H5parm and beam
      const MC2x2F gain_combined_1 = MultiplyGains<Mode>(gain_h5_1, gain_b_1);
      const MC2x2F gain_combined_2 = MultiplyGains<Mode>(gain_h5_2, gain_b_2);

      if (apply_forward) {
        ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
      }
      ApplyConjugatedGain<Mode>(data, gain_combined_1, gain_combined_2);

      const float weighted_squared_gain_h5 =
          ComputeWeightedSquaredGain<Mode>(gain_h5_1, gain_h5_2, weights);
      h5Sum += weighted_squared_gain_h5 * image_weights[ch];

      const float weighted_squared_gain_combined =
          ComputeWeightedSquaredGain<Mode>(gain_combined_1, gain_combined_2,
                                           weights);

      correctionSum += weighted_squared_gain_combined * image_weights[ch];

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Apply facet beam
      const size_t beam_offset = ch * _pointResponseBufferSize;
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);
      if (apply_forward) {
        ApplyGain<Mode>(data, gain_b_1, gain_b_2);
      }
      ApplyConjugatedGain<Mode>(data, gain_b_1, gain_b_2);
      // Apply h5
      // Column major indexing
      const size_t offset_h5 =
          (_timeOffset[ms_index] * n_channels + ch) * n_stations * nparms;
      const size_t offset_h5_1 = offset_h5 + antenna1 * nparms;
      const size_t offset_h5_2 = offset_h5 + antenna2 * nparms;
      const MC2x2F gain_h5_1(&_cachedParmResponse[ms_index][offset_h5_1]);
      const MC2x2F gain_h5_2(&_cachedParmResponse[ms_index][offset_h5_2]);
      if (apply_forward) {
        ApplyGain<Mode>(data, gain_h5_1, gain_h5_2);
      }
      ApplyConjugatedGain<Mode>(data, gain_h5_1, gain_h5_2);
      const float weighted_squared_gain_h5 =
          ComputeWeightedSquaredGain<Mode>(gain_h5_1, gain_h5_2, weights);
      h5Sum += weighted_squared_gain_h5 * image_weights[ch];

      const MC2x2F gain_combined_1(gain_b_1[0] * gain_h5_1[0], 0, 0,
                                   gain_b_1[3] * gain_h5_1[3]);
      const MC2x2F gain_combined_2(gain_b_2[0] * gain_h5_2[0], 0, 0,
                                   gain_b_2[3] * gain_h5_2[3]);

      const float weighted_squared_gain_combined =
          ComputeWeightedSquaredGain<Mode>(gain_combined_1, gain_combined_2,
                                           weights);

      correctionSum += weighted_squared_gain_combined * image_weights[ch];

      data += GetNVisibilities(Mode);
      weights += GetNVisibilities(Mode);
    }
  }
  correction_sum_ += correctionSum;
  h5_correction_sum_ += h5Sum;
}

template void VisibilityModifier::ApplyConjugatedDual<GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<GainMode::kTrace>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<GainMode::k2VisDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

template void VisibilityModifier::ApplyConjugatedDual<GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward);

#endif  // HAVE_EVERYBEAM
