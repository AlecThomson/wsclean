#ifndef RENDERER_H_
#define RENDERER_H_

#include <cstring>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include "../model/model.h"

namespace renderer {

/**
 * @brief Struct collecting the relevant image coordinate settings for
 * rendereing the source components.
 *
 */
struct ImageCoordinateSettings {
  ImageCoordinateSettings() = default;

  /**
   * @brief Extract coordinate settings from an aocommon::FitsReader object.
   *
   * @param fits_reader aocommon::FitsReader object.
   */
  ImageCoordinateSettings(const aocommon::FitsReader& fits_reader)
      : phase_centre_ra(fits_reader.PhaseCentreRA()),
        phase_centre_dec(fits_reader.PhaseCentreDec()),
        pixel_scale_l(fits_reader.PixelSizeX()),
        pixel_scale_m(fits_reader.PixelSizeY()),
        phase_centre_dl(fits_reader.PhaseCentreDL()),
        phase_centre_dm(fits_reader.PhaseCentreDM()) {}

  long double phase_centre_ra;
  long double phase_centre_dec;
  long double pixel_scale_l;
  long double pixel_scale_m;
  long double phase_centre_dl;
  long double phase_centre_dm;
};

/**
 * Restore model component using a circular beam, goverened by
 * \c beam_size .
 */
void RestoreWithCircularBeam(aocommon::Image& image,
                             const ImageCoordinateSettings& image_settings,
                             const Model& model, long double beam_size,
                             long double start_frequency,
                             long double end_frequency,
                             aocommon::PolarizationEnum polarization);

/**
 * Restore a model image with an elliptical beam
 */
void RestoreWithEllipticalBeam(
    aocommon::Image& image, const ImageCoordinateSettings& image_settings,
    const Model& model, long double beam_major_axis,
    long double beam_minor_axis, long double beam_position_angle,
    long double start_frequency, long double end_frequency,
    aocommon::PolarizationEnum polarization, size_t thread_count);

// class Renderer {
//  public:
//   Renderer(const Model& model, long double phase_centre_ra,
//            long double phase_centre_dec, long double pixel_scale_l,
//            long double pixel_scale_m, long double phase_centre_dl,
//            long double phase_centre_dm)
//       : model_(model),
//         phase_centre_ra_(phase_centre_ra),
//         phase_centre_dec_(phase_centre_dec),
//         pixel_scale_l_(pixel_scale_l),
//         pixel_scale_m_(pixel_scale_m),
//         phase_centre_dl_(phase_centre_dl),
//         phase_centre_dm_(phase_centre_dm) {}

//   Renderer(const Renderer&) = delete;
//   Renderer& operator=(const Renderer&) = delete;

//   /**
//    * Restore model component using a circular beam, goverened by
//    * \c beam_size .
//    */
//   void RestoreWithCircularBeam(double* image_data, size_t image_width,
//                                size_t image_height, const Model& model,
//                                long double beam_size,
//                                long double start_frequency,
//                                long double end_frequency,
//                                aocommon::PolarizationEnum polarization)
//                                const;

//   /**
//    * Restore a model with an elliptical beam
//    */
//   void RestoreWithEllipticalBeam(float* image_data, size_t image_width,
//                                  size_t image_height, long double beamMaj,
//                                  long double beamMin, long double beamPA,
//                                  long double start_frequency,
//                                  long double end_frequency,
//                                  aocommon::PolarizationEnum polarization,
//                                  size_t thread_count) const;

//  private:
//   /**
//    * Render without beam convolution, such that each point-source is one
//    pixel.
//    */
//   void RenderModel(float* image_data, size_t image_width, size_t
//   image_height,
//                    long double start_frequency, long double end_frequency,
//                    aocommon::PolarizationEnum polarization) const;

//   const Model& model_;
//   const long double phase_centre_ra_;
//   const long double phase_centre_dec_;
//   const long double pixel_scale_l_;
//   const long double pixel_scale_m_;
//   const long double phase_centre_dl_;
//   const long double phase_centre_dm_;
// };
}  // namespace renderer

#endif
