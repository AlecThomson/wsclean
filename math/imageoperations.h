#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include "../structures/outputchannelinfo.h"
#include "../io/imagefilename.h"

#include <aocommon/hmatrix4x4.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <vector>
#include <string>

class Settings;

namespace math {

/**
 * Multiplies every pixel of a set of Stokes I,Q,U,V images with a Mueller
 * matrix. This function assumes that the Mueller matrix is a correction to
 * linear instrumental polarizations, i.e. it applies a correction on
 * XX/XY/YX/YY values. The Stokes IQUV image values are therefore first
 * converted to linear polarization, then the matrix is applied, and then the
 * values are converted and written back to Stokes images.
 */
void CorrectImagesForMuellerMatrix(const aocommon::HMC4x4& mueller_correction,
                                   std::array<aocommon::Image*, 4>& images);

}  // namespace math

class ImageOperations {
 public:
  static void FitBeamSize(const Settings& settings, double& bMaj, double& bMin,
                          double& bPA, const aocommon::Image& image,
                          double beamEstimate);

  static void DetermineBeamSize(const Settings& settings, double& bMaj,
                                double& bMin, double& bPA, double& bTheoretical,
                                const aocommon::Image& image,
                                double initialEstimate);

  static void MakeMFSImage(const Settings& settings,
                           const std::vector<OutputChannelInfo>& infoPerChannel,
                           OutputChannelInfo& mfsInfo,
                           const std::string& suffix, size_t intervalIndex,
                           aocommon::PolarizationEnum pol,
                           ImageFilenameType image_type,
                           std::optional<size_t> directionIndex = std::nullopt);

  static void RenderMFSImage(const Settings& settings,
                             const OutputChannelInfo& mfsInfo,
                             size_t intervalIndex,
                             aocommon::PolarizationEnum pol, bool isImaginary,
                             bool isPBCorrected);
};

#endif
