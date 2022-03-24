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
 * @brief Restore a model image with an elliptical beam.
 *
 * @param image Image to which restored sources are written.
 * @param image_settings Image coordinate settings.
 * @param model Modeled source components.
 * @param beam_major_axis Major axis of elliptical beam to be applied [rad].
 * @param beam_minor_axis Minor axis of elliptical beam to be applied [rad].
 * @param beam_position_angle Position angle of beam [rad].
 * @param start_frequency Start frequency [Hz].
 * @param end_frequency End frequency [Hz].
 * @param polarization Polarization enum.
 * @param thread_count Numbers of threads to use.
 */
void RestoreWithEllipticalBeam(
    aocommon::Image& image, const ImageCoordinateSettings& image_settings,
    const Model& model, long double beam_major_axis,
    long double beam_minor_axis, long double beam_position_angle,
    long double start_frequency, long double end_frequency,
    aocommon::PolarizationEnum polarization, size_t thread_count);

}  // namespace renderer

#endif
