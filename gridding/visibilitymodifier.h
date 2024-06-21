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

#include "averagecorrection.h"
#include "gainmode.h"

class SynchronizedMS;

/**
 * Control whether to apply the modifier, sum the corrections or both when
 * calling methods that apply a modifier e.g. @ref ApplyConjugatedParmResponse()
 */
enum class ModifierBehaviour { kApply, kSum, kApplyAndSum };

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

  void ResetCache(size_t n_measurement_sets) {
    // Assign, rather than a resize here to make sure that
    // caches are re-initialized - even in the case an MSGridderBase
    // object would be re-used for multiple gridding tasks.
    _cachedParmResponse.clear();
    _cachedMSTimes.clear();
    _timeOffsets.clear();
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
   * @tparam Behaviour Determines whether we should apply the gain, sum the
   * correction or both.
   * @tparam Mode Gain application mode that defines how the gain is
   * applied, the PolarizationCount is also implied/determined by the gain mode.
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
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedParmResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights, size_t ms_index,
                                   size_t n_channels, size_t n_antennas,
                                   size_t antenna1, size_t antenna2,
                                   bool apply_forward);

  template <GainMode GainEntry>
  void ApplyParmResponse(std::complex<float>* data, size_t ms_index,
                         size_t n_channels, size_t n_antennas, size_t antenna1,
                         size_t antenna2);

  void SetMSTimes(size_t ms_index, std::vector<double>&& times) {
    _cachedMSTimes[ms_index] = std::move(times);
  }

  void SetH5Parm(
      const std::vector<schaapcommon::h5parm::H5Parm>& h5parms,
      const std::vector<schaapcommon::h5parm::SolTab*>& first_solutions,
      const std::vector<schaapcommon::h5parm::SolTab*>& second_solutions,
      const std::vector<schaapcommon::h5parm::GainType>& gain_types) {
    _h5parms = &h5parms;
    _firstSolutions = &first_solutions;
    _secondSolutions = &second_solutions;
    _gainTypes = &gain_types;
  }

  bool HasH5Parm() const { return _h5parms && !_h5parms->empty(); }

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Compute and cache the beam response if no cached response
   * present for the provided time.
   */
  void CacheBeamResponse(double time, size_t fieldId,
                         const aocommon::BandData& band);

  template <GainMode Mode>
  void ApplyBeamResponse(std::complex<float>* data, size_t n_channels,
                         size_t antenna1, size_t antenna2);

  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedBeamResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights,
                                   size_t n_channels, size_t antenna1,
                                   size_t antenna2, bool apply_forward);

  /**
   * Correct the data for both the conjugated beam and the
   * conjugated h5parm solutions.
   */
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedDual(std::complex<float>* data, const float* weights,
                           const float* image_weights, size_t n_channels,
                           size_t n_stations, size_t antenna1, size_t antenna2,
                           size_t ms_index, bool apply_forward);
#endif

  void SetFacetDirection(double ra, double dec) {
    _facetDirectionRA = ra;
    _facetDirectionDec = dec;
  }
  double FacetDirectionRA() const { return _facetDirectionRA; }
  double FacetDirectionDec() const { return _facetDirectionDec; }
  void ResetSums() {
    correction_sum_ = AverageCorrection();
    beam_correction_sum_ = AverageCorrection();
  }
  /**
   * Sum of the full corrections applied to the visibilities. In case both
   * beam and h5parm solutions are applied, this is the weighed sum over the
   * (squared) product of both. Otherwise it is over the (squared) contribution
   * of either the beam or solutions.
   *
   * It is not an average yet, because this class doesn't store the sum of
   * weights. It would be a redundant calculation, because the gridder
   * already calculates the sum of weights.
   */
  const AverageCorrection& TotalCorrectionSum() const {
    return correction_sum_;
  }
  /**
   * In case both beam and solution gains are applied, this represents the beam
   * part of the corrections. In that case, it is the weighed sum of the squared
   * beam matrices. This is used in the final primary beam correction to
   * (approximately) separate the beam and solution parts to be able to apply a
   * smooth beam correction. If only one correction is applied, this value
   * remains zero.
   *
   * Like @ref TotalCorrectionSum(), this should be divided by the sum of
   * weights to turn it into the average correction.
   */
  const AverageCorrection& BeamCorrectionSum() const {
    return beam_correction_sum_;
  }

  /**
   * Number of complex values per solution: 4 in the case of fulljones,
   * otherwise 2 (scalar solutions are already converted to dual solutions
   * by the h5parm reader).
   */
  inline size_t NValuesPerSolution(size_t ms_index) const {
    using schaapcommon::h5parm::GainType;
    const size_t solution_index = (*_gainTypes).size() == 1 ? 0 : ms_index;
    return (*_gainTypes)[solution_index] == GainType::kFullJones ? 4 : 2;
  }

 private:
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
   * Element ms_index is a vector of complex gains of
   * size n_times x n_channels x n_stations x n_parameters, where n_parameters
   * is the fastest changing. n_parameters is 2 (for diagonal) or
   * 4 (for full jones).
   *
   * The ms_index may not be a consecutive index because this gridder
   * may not have to grid all measurement sets that were specified to
   * wsclean.
   */
  std::map<size_t, std::vector<std::complex<float>>> _cachedParmResponse;
  /**
   * Each element holds a vector with the measurement set times. The map
   * is indexed by a (non-consecutive) ms_index.
   */
  std::map<size_t, std::vector<double>> _cachedMSTimes;
  /**
   * Each element holds the current offset position into the _cachedParmResponse
   * and _cachedMSTimes elements of the same ms_index. The map is indexed by a
   * (non-consecutive) ms_index.
   */
  std::map<size_t, size_t> _timeOffsets;
  /**
   * Optional pointers to vectors with h5parm solution objects.
   * The GriddingTaskManager, which always outlives GriddingTasks and their
   * VisibilityModifier, owns the objects.
   * If all measurement sets use the same solution, the vectors have one
   * element. Otherwise, they have one element for each ms.
   * The second solution table is optional. If both tables are used, the first
   * table has the amplitude part and the second table has the phase part.
   * @{
   */
  const std::vector<schaapcommon::h5parm::H5Parm>* _h5parms = nullptr;
  const std::vector<schaapcommon::h5parm::SolTab*>* _firstSolutions = nullptr;
  const std::vector<schaapcommon::h5parm::SolTab*>* _secondSolutions = nullptr;
  const std::vector<schaapcommon::h5parm::GainType>* _gainTypes = nullptr;
  /** @} */
  size_t _pointResponseBufferSize = 0;
  double _facetDirectionRA = 0.0;
  double _facetDirectionDec = 0.0;
  AverageCorrection correction_sum_;
  AverageCorrection beam_correction_sum_;
};

#endif
