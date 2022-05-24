#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include "../structures/outputchannelinfo.h"

#include <aocommon/polarization.h>

#include <vector>
#include <string>

class ImageOperations {
 public:
  static void FitBeamSize(const class Settings& settings, double& bMaj,
                          double& bMin, double& bPA, const float* image,
                          double beamEstimate);

  static void DetermineBeamSize(const class Settings& settings, double& bMaj,
                                double& bMin, double& bPA, double& bTheoretical,
                                const float* image, double initialEstimate);

  static void MakeMFSImage(const class Settings& settings,
                           const std::vector<OutputChannelInfo>& infoPerChannel,
                           OutputChannelInfo& mfsInfo,
                           const std::string& suffix, size_t intervalIndex,
                           aocommon::PolarizationEnum pol, bool isImaginary,
                           bool isPSF = false);

  static void RenderMFSImage(const class Settings& settings,
                             const OutputChannelInfo& mfsInfo,
                             size_t intervalIndex,
                             aocommon::PolarizationEnum pol, bool isImaginary,
                             bool isPBCorrected);

  static void RemoveNansInPsf(float* psf, size_t width, size_t height) {
    float* end = psf + width * height;
    while (psf != end) {
      if (!std::isfinite(*psf)) *psf = 0.0;
      ++psf;
    }
  }
};

#endif
