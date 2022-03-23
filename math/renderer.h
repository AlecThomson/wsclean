#ifndef RENDERER_H_
#define RENDERER_H_

#include <cstring>

#include <aocommon/polarization.h>

#include "../model/model.h"

class Renderer {
 public:
  Renderer(const Model& model, long double phase_centre_ra,
           long double phase_centre_dec, long double pixel_scale_l,
           long double pixel_scale_m, long double phase_centre_dl,
           long double phase_centre_dm)
      : model_(model),
        phase_centre_ra_(phase_centre_ra),
        phase_centre_dec_(phase_centre_dec),
        pixel_scale_l_(pixel_scale_l),
        pixel_scale_m_(pixel_scale_m),
        phase_centre_dl_(phase_centre_dl),
        phase_centre_dm_(phase_centre_dm) {}

  Renderer(const Renderer&) = delete;
  Renderer& operator=(const Renderer&) = delete;

  /**
   * Restore model component using a circular beam, goverened by
   * \c beamSize .
   */
  void RestoreWithCircularBeam(double* image_data, size_t image_width,
                               size_t image_height, const Model& model,
                               long double beamSize,
                               long double start_frequency,
                               long double end_frequency,
                               aocommon::PolarizationEnum polarization) const;

  /**
   * Restore a model with an elliptical beam
   */
  void RestoreWithEllipticalBeam(float* image_data, size_t image_width,
                                 size_t image_height, long double beamMaj,
                                 long double beamMin, long double beamPA,
                                 long double start_frequency,
                                 long double end_frequency,
                                 aocommon::PolarizationEnum polarization,
                                 size_t thread_count) const;

 private:
  /**
   * Render without beam convolution, such that each point-source is one
   pixel.
   */
  void RenderModel(float* image_data, size_t image_width, size_t image_height,
                   long double start_frequency, long double end_frequency,
                   aocommon::PolarizationEnum polarization) const;

  const Model& model_;
  const long double phase_centre_ra_;
  const long double phase_centre_dec_;
  const long double pixel_scale_l_;
  const long double pixel_scale_m_;
  const long double phase_centre_dl_;
  const long double phase_centre_dm_;
};

#endif
