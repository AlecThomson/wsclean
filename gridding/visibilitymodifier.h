#ifndef GRIDDING_VISIBILITY_MODIFIER_H_
#define GRIDDING_VISIBILITY_MODIFIER_H_

#include <stdexcept>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <schaapcommon/h5parm/jonesparameters.h>

#include <aocommon/banddata.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

class SynchronizedMS;

/**
 * Enum for selecting the entry or entries from the direction dependent gain
 * matrix that are to be used to correct the visibilities during the reading
 * and/or writing operations.
 */
enum class GainMode {
  /// Correct visibilities only with the X solution.
  kXX,
  /// Correct visibilities only with the Y solution.
  kYY,
  /// Correct visibilities with the X and Y solutions, but not with the
  /// cross-terms (XY/YX).
  /// When the input is Stokes I visibilities, the trace of the gain is taken
  /// and applied. When the input contains both XX and YY or all four
  /// polarizations, the X/Y gains are separately applied.
  kDiagonal,
  /// Correct visibilities with the full 2x2 complex matrix.
  kFull
};

inline GainMode GetGainMode(aocommon::PolarizationEnum polarization,
                            size_t n_visibility_polarizations) {
  switch (n_visibility_polarizations) {
    case 1:
      switch (polarization) {
        case aocommon::Polarization::XX:
          return GainMode::kXX;
        case aocommon::Polarization::YY:
          return GainMode::kYY;
        case aocommon::Polarization::StokesI:
          // polarization might also be RR. We still need to provide a GainMode,
          // so we also return diagonal in those cases.
        default:
          return GainMode::kDiagonal;
      }
      break;
    case 2:
      if (polarization == aocommon::Polarization::DiagonalInstrumental ||
          polarization == aocommon::Polarization::StokesI)
        return GainMode::kDiagonal;
      break;
    case 4:
      if (polarization == aocommon::Polarization::DiagonalInstrumental)
        return GainMode::kDiagonal;
      else if (polarization == aocommon::Polarization::Instrumental)
        return GainMode::kFull;
      break;
  }
  throw std::runtime_error(
      "Invalid combination of polarization (" +
      aocommon::Polarization::TypeToFullString(polarization) +
      ") and n_visibility_polarizations (" +
      std::to_string(n_visibility_polarizations) + ") in GetGainMode()");
}

/**
 * Applies beam and h5parm solutions to visibilities.
 * See the documentation for function @ref ApplyConjugatedParmResponse()
 * for an overview of parameters that hold for most of these functions.
 */
class VisibilityModifier {
 public:
  VisibilityModifier() = default;

  void InitializePointResponse(SynchronizedMS&& ms,
                               double facet_beam_update_time,
                               const std::string& element_response,
                               size_t n_channels,
                               const std::string& data_column,
                               const std::string& mwa_path);

  /**
   * A function that initializes this visibility modifier for testing. After
   * calling this function, the Apply* functions can be called.
   * @param beam_response vector of n_channels * n_stations * 4 gain elements.
   * @param parm_response vector of n_channels * n_stations * 2
   */
  void InitializeMockResponse(
      size_t n_channels, size_t n_stations,
      const std::vector<std::complex<double>>& beam_response,
      const std::vector<std::complex<float>>& parm_response);

  void SetNoPointResponse() {
#ifdef HAVE_EVERYBEAM
    _pointResponse = nullptr;
    _cachedBeamResponse.clear();
#endif
  }

  void SetBeamInfo(std::string mode, std::string normalisation) {
#ifdef HAVE_EVERYBEAM
    _beamModeString = std::move(mode);
    _beamNormalisationMode = std::move(normalisation);
#endif
  }

  void ResetCache(size_t n_measurement_sets,
                  const std::vector<std::string>& solutionFiles,
                  const std::vector<std::string>& solutionTables) {
    // Assign, rather than a resize here to make sure that
    // caches are re-initialized - even in the case an MSGridderBase
    // object would be re-used for multiple gridding tasks.
    _cachedParmResponse.assign(n_measurement_sets, {});
    _cachedMSTimes.assign(n_measurement_sets, {});
    _timeOffset.assign(n_measurement_sets, 0u);
    SetH5Parms(n_measurement_sets, solutionFiles, solutionTables);
  }

  /**
   * @brief Cache the solutions from a h5 solution file and update the
   * associated time.
   */
  void CacheParmResponse(double time,
                         const std::vector<std::string>& antennaNames,
                         const aocommon::BandData& band, size_t ms_index);

  /**
   * Applies the conjugate (is backward, or imaging direction) h5parm gain
   * to given data.
   * @tparam PolarizationCount Number of differently-polarized correlations in
   * the data, e.g. 2 for XX/YY and 4 for full Jones matrices.
   * @tparam GainEntry Gain application mode that defines how the gain is
   * applied.
   * @param [in,out] data Data array with n_channels x PolarizationCount
   * elements.
   * @param weights Array with for each data value the corresponding weight.
   * @param image_weights Array of size n_channels (polarizations are assumed to
   * have equal imaging weights) with the imaging weighting mode (e.g. from
   * Briggs weighting).
   * @param apply_forward If true, an additional (non-conjugate) forward gain is
   * applied to the data. This is necessary for calculating direction-dependent
   * PSFs.
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  float ApplyConjugatedParmResponse(std::complex<float>* data,
                                    const float* weights,
                                    const float* image_weights, size_t ms_index,
                                    size_t n_channels, size_t n_antennas,
                                    size_t antenna1, size_t antenna2,
                                    bool apply_forward);

  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyParmResponse(std::complex<float>* data, size_t ms_index,
                         size_t n_channels, size_t n_antennas, size_t antenna1,
                         size_t antenna2);

  void SetMSTimes(size_t ms_index, std::vector<double>&& times) {
    _cachedMSTimes[ms_index] = std::move(times);
  }

  bool HasH5Parm() const { return !_h5parms.empty(); }

  struct DualResult {
    float h5Sum = 0.0;
    float correctionSum = 0.0;
  };

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Compute and cache the beam response if no cached response
   * present for the provided time.
   */
  void CacheBeamResponse(double time, size_t fieldId,
                         const aocommon::BandData& band);

  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyBeamResponse(std::complex<float>* data, size_t n_channels,
                         size_t antenna1, size_t antenna2);

  template <size_t PolarizationCount, GainMode GainEntry>
  float ApplyConjugatedBeamResponse(std::complex<float>* data,
                                    const float* weights,
                                    const float* image_weights,
                                    size_t n_channels, size_t antenna1,
                                    size_t antenna2, bool apply_forward);

  /**
   * Correct the data for both the conjugated beam and the
   * conjugated h5parm solutions.
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  DualResult ApplyConjugatedDual(std::complex<float>* data,
                                 const float* weights,
                                 const float* image_weights, size_t n_channels,
                                 size_t n_stations, size_t antenna1,
                                 size_t antenna2, size_t ms_index,
                                 bool apply_forward);
#endif

  void SetFacetDirection(double ra, double dec) {
    _facetDirectionRA = ra;
    _facetDirectionDec = dec;
  }
  double FacetDirectionRA() const { return _facetDirectionRA; }
  double FacetDirectionDec() const { return _facetDirectionDec; }

 private:
  void SetH5Parms(size_t n_measurement_sets,
                  const std::vector<std::string>& solutionFiles,
                  const std::vector<std::string>& solutionTables);

#ifdef HAVE_EVERYBEAM
  // _telescope attribute needed to keep the telecope in _pointResponse alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  /**
   * The beam response for the currently processed timestep.
   * It's of size n_channels x _pointResponseBufferSize, which equals
   * n_channels x n_stations x n_elements(=4), where n_elements is the fastest
   * changing index.
   */
  aocommon::UVector<std::complex<float>> _cachedBeamResponse;
  everybeam::BeamMode _beamMode = everybeam::BeamMode::kNone;
#endif
  std::string _beamModeString;
  std::string _beamNormalisationMode;
  /**
   * _cachedParmResponse[ms_index] is a vector of complex gains of
   * size n_times x n_channels x n_stations x n_parameters, where n_parameters
   * is the fastest changing. n_parameters is 2 (for diagonal) or
   * 4 (for full jones).
   */
  std::vector<std::vector<std::complex<float>>> _cachedParmResponse;
  std::vector<std::unique_ptr<schaapcommon::h5parm::H5Parm>> _h5parms;
  std::vector<
      std::pair<schaapcommon::h5parm::SolTab*, schaapcommon::h5parm::SolTab*>>
      _h5SolTabs;
  /// The correction type for each measurement set (given ms_index)
  std::vector<schaapcommon::h5parm::JonesParameters::CorrectType>
      _h5correctType;
  std::vector<std::vector<double>> _cachedMSTimes;
  std::vector<size_t> _timeOffset;
  size_t _pointResponseBufferSize = 0;
  double _facetDirectionRA = 0.0;
  double _facetDirectionDec = 0.0;
};

#endif
