#ifndef WSCLEAN_DECONVOLUTION_SETTINGS_H_
#define WSCLEAN_DECONVOLUTION_SETTINGS_H_

#include <set>

#include <aocommon/polarization.h>
#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "../multiscale/multiscaletransforms.h"

/**
 * @brief The value of LocalRmsMethod describes how the RMS map
 * should be used.
 */
enum class LocalRmsMethod { kNone, kRmsWindow, kRmsAndMinimumWindow };

struct DeconvolutionSettings {
  enum class Algorithm { kPython, kMoreSane, kIuwt, kMultiScale, kGeneric };

  /**
   * @{
   * Settings that are duplicates from top level settings, and also used outside
   * deconvolution.
   */
  size_t trimmedImageWidth = 0;
  size_t trimmedImageHeight = 0;
  size_t channelsOut = 1;
  double pixelScaleX = 0.0;
  double pixelScaleY = 0.0;
  size_t threadCount = aocommon::system::ProcessorCount();
  std::string prefixName = "wsclean";
  /** @} */

  /**
   * @{
   * These settings strictly pertain to deconvolution only.
   */
  std::set<aocommon::PolarizationEnum> linkedPolarizations;
  size_t parallelDeconvolutionMaxSize = 0;
  size_t parallelDeconvolutionMaxThreads = 0;
  double deconvolutionThreshold = 0.0;
  double deconvolutionGain = 0.1;
  double deconvolutionMGain = 1.0;
  bool autoDeconvolutionThreshold = false;
  bool autoMask = false;
  double autoDeconvolutionThresholdSigma = 0.0;
  double autoMaskSigma = 0.0;
  LocalRmsMethod localRMSMethod = LocalRmsMethod::kNone;
  double localRMSWindow = 25.0;
  std::string localRMSImage;
  bool saveSourceList = false;
  size_t deconvolutionIterationCount = 0;
  size_t majorIterationCount = 0;
  bool allowNegativeComponents = true;
  bool stopOnNegativeComponents = false;
  bool squaredJoins = false;
  double spectralCorrectionFrequency = 0.0;
  std::vector<float> spectralCorrection;
  double deconvolutionBorderRatio = 0.0;
  std::string fitsDeconvolutionMask;
  std::string casaDeconvolutionMask;
  bool horizonMask = false;
  double horizonMaskDistance = 0.0;
  schaapcommon::fitters::SpectralFittingMode spectralFittingMode =
      schaapcommon::fitters::SpectralFittingMode::NoFitting;
  size_t spectralFittingTerms = 0;
  std::string forcedSpectrumFilename;
  /**
   * The number of channels used during deconvolution. This can be used to
   * image with more channels than deconvolution. Before deconvolution,
   * channels are averaged, and after deconvolution they are interpolated.
   * If it is 0, all channels should be used.
   */
  size_t deconvolutionChannelCount = 0;

  Algorithm algorithm = Algorithm::kGeneric;

  struct {
    std::string deconvolutionFilename;
  } python;

  struct {
    std::string location;
    std::string args;
    std::vector<double> sigmaLevels;
  } moreSane;

  struct {
    bool SNRTest = false;
  } iuwt;

  struct {
    bool fastSubMinorLoop = true;
    double gain = 0.2;
    double deconvolutionScaleBias = 0.6;
    size_t maxScales = 0;
    double convolutionPadding = 1.1;
    std::vector<double> scaleList;
    MultiScaleTransforms::Shape shapeFunction =
        MultiScaleTransforms::TaperedQuadraticShape;
  } multiscale;

  struct {
    bool useSubMinorOptimization = true;
  } generic;
  /** @} */
};

#endif
